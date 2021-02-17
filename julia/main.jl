include("logger.jl")
include("helper/p3.jl")
include("helper/graph.jl")
include("heuristik/main.jl")
#include("gurobi/main.jl")
include("gurobi/pythonWrapper.jl")
include("data_reduction/data_reduction.jl")
include("data_reduction/merge_greater_k.jl")
include("data_reduction/heavy_edge.jl")
include("data_reduction/neighborhood.jl")

P3List = Array{P3,1}

recursive_steps = 0

function applyDataReductionBranch(g::Graph, c::Float64, ub::Float64, params::Params, recDepth::Int)
    k_decr = 0.
    path::Array{Operation, 1} = []

    he_be = params.applyHeavyBothEnds && recDepth % params.applyHeavyBothEndsFreq == 0
    he_se = params.applyHeavySingleEnd && recDepth % params.applyHeavySingleEndFreq == 0
    hne = params.applyHeavyNonEdge && recDepth % params.applyHeavyNonEdgeFreq == 0

    if he_be || he_se || hne
        be_path, be_k_decr = heavyEdgeEndsRec(g, he_be, he_se, hne)
        k_decr += be_k_decr
        path = vcat(path, be_path)
    end

    if params.applyMergeGreaterK && recDepth % params.applyMergeGreaterKFreq == 0
        mg_path, mg_k_decr = applyMergeGreaterKRec(g, ub - c)
        k_decr += mg_k_decr
        path = vcat(path, mg_path)
    end

    if params.applyLN && recDepth % params.applyLNFreq == 0
        ln_path, ln_k_decr = large_neighborhood_rec(g)
        k_decr += ln_k_decr
        path = vcat(path, ln_path)
    end

    reverse!(path)

    return path, k_decr
end

function solve(g::Graph, c::Float64, ub::Float64, params::Params, noP3s::P3List=[], existingP3s::P3List=[], start::Int=1, recDepth::Int=0)::Tuple{Float64, Array{Operation, 1}}
    global recursive_steps += 1

    dataReductionPath, k_decr = applyDataReductionBranch(g, c, ub, params, recDepth)
    c += k_decr

    if params.applyLB2 && recDepth % params.applyLB2Freq == 0
        lb2 = getMinWeightsOfNOP3s(g, noP3s)
        if c + lb2 > ub
            undoK(g, size(dataReductionPath, 1))
            return ub, dataReductionPath
        end
    end

    if params.applyLB3 && recDepth % params.applyLB3Freq == 0
        lb3 = gurobiLB(g)
        if c + lb3 > ub
            undoK(g, size(dataReductionPath, 1))
            return ub, dataReductionPath
        end
    end

    if params.applyComponentSplitting && recDepth % params.applyComponentSplittingFreq == 0 && !isConnected(g)
        comps = getComponentGraphs(g)
        path = []
        for comp in comps
            # IDEE: UB neu berechnen für comp?
            c_temp, path_temp = solve(comp, 0., ub - c, params, noP3s, existingP3s, 1, recDepth)
            c += c_temp
            path = vcat(path, path_temp)
        end

        if ub < c
            undoK(g, size(dataReductionPath, 1))
            return ub, dataReductionPath
        else
            undoK(g, size(dataReductionPath, 1))
            path = vcat(path, dataReductionPath)
            return c, path
        end
    end

    foundP3 = false
    if start > 0
        s = size(existingP3s, 1)
        for edge_idx = start:s
            @inbounds p = existingP3s[edge_idx]
            if @inbounds g.indexes_lookup[p[1]] &&
               @inbounds g.indexes_lookup[p[2]] &&
               @inbounds g.indexes_lookup[p[3]] &&
               checkP3(g, p)
                foundP3 = true
                start = edge_idx
                break
            end
        end
    end

    if !foundP3
        foundP3, existingP3s = findP3AdjList(g)

        noP3s = params.applyLB2 ? findNonOverlappingP3s(g, existingP3s) : P3[]
        start = 1
    end

    if !foundP3
        undoK(g, size(dataReductionPath, 1))
        return c, dataReductionPath
    else
        s, m, e = existingP3s[start]

        a = s
        b = m
        if getWeight(g, m, e) > getWeight(g, s, m)
            a = m
            b = e
        end

        deleteCost::Int = -1
        deletePath::Array{Operation, 1} = []
        mergeCost::Int = -1
        mergePath::Array{Operation, 1} = []

        if params.applyBranchDeleteFirst
            # DELETE ===========
            weight = getWeight(g, a, b)

            setWeight(g, a, b, -Inf)
            deleteCost, deletePath = solve(g, c + abs(weight), ub, params, noP3s, existingP3s, start + 1, recDepth + 1)
            setWeight(g, a, b, weight)

            ub = min(ub, deleteCost)
        end

        # MERGE ============

        k_decr = merge!(g, a, b)
        mergeCost, mergePath = solve(g, c + k_decr, ub, params, noP3s, existingP3s, start + 1, recDepth + 1)
        undoK(g, 1)

        if !params.applyBranchDeleteFirst
            ub = min(ub, mergeCost)

            # DELETE ===========
            weight = getWeight(g, a, b)

            setWeight(g, a, b, -Inf)
            deleteCost, deletePath = solve(g, c + abs(weight), ub, params, noP3s, existingP3s, start + 1, recDepth + 1)
            setWeight(g, a, b, weight)
        end

        undoK(g, size(dataReductionPath, 1))

        if params.applyBranchDeleteFirst
            if deleteCost <= mergeCost
                push!(deletePath, Operation(deleteOp, a, b))
                deletePath = vcat(deletePath, dataReductionPath)
                return deleteCost, deletePath
            else
                push!(mergePath, Operation(mergeOp, a, b))
                mergePath = vcat(mergePath, dataReductionPath)
                return mergeCost, mergePath
            end
        else
            if deleteCost < mergeCost
                push!(deletePath, Operation(deleteOp, a, b))
                deletePath = vcat(deletePath, dataReductionPath)
                return deleteCost, deletePath
            else
                push!(mergePath, Operation(mergeOp, a, b))
                mergePath = vcat(mergePath, dataReductionPath)
                return mergeCost, mergePath
            end
        end
    end
end

function getMinWeightsOfNOP3s(g::Graph, p3s::P3List)::Int
    sum = 0.
    for p = p3s
        if @inbounds g.indexes_lookup[p[1]] &&
            @inbounds g.indexes_lookup[p[2]] &&
            @inbounds g.indexes_lookup[p[3]] &&
            checkP3(g, p)
            sum += min(getWeight(g, p[1], p[2]), getWeight(g, p[2], p[3]), -getWeight(g, p[1], p[3]))
        end
    end
    return sum
end

function findNonOverlappingP3s(g::Graph, p3s::P3List)::P3List
    marked = fill(false, (g.n_total, g.n_total))

    ret = []
    for p3 in p3s
        s, m, e = p3
        if !marked[s, m] && !marked[m, e] && !marked[s, e]
            push!(ret, p3)
            marked[s, m] = true
            marked[m, s] = true
            marked[m, e] = true
            marked[e, m] = true
            marked[s, e] = true
            marked[e, s] = true
        end
    end

    return ret
end

function solveInstance(g::Graph, params::Params)
    g_before = copy(g)

    k_decr::Float64 = doDataReductionPre(g)

    println("#k_decr: $k_decr")

    # If there are less than 3 nodes it is a cluster graph
    if size(g.indexes, 1) < 3
        g.foundSolution = true
        undoAll!(g)
        printSolution(g_before, g)
        println("#Finished: DR")
        return k_decr
    end

    lowerBound = gurobiLB(g) + k_decr
    upperBoundGraph, upperBound = solveHeuristik(g_before, params) # Check with DR

    println("#Initial LB: $lowerBound UB: $upperBound")

    # If lowerBound == upperBound we found an optimal solution
    if lowerBound == upperBound
        printSolution(g_before, upperBoundGraph)
        println("#Finished: eq")
        return upperBound
    end

    # Search
    found, p3s = findP3AdjList(g)

    nonOverlappingP3s = params.applyLB2Freq != 0 ? findNonOverlappingP3s(g, p3s) : P3[]

    minCost, path = solve(g, k_decr, upperBound, params, nonOverlappingP3s, p3s)

    if minCost == upperBound
        #printEdges(upperBoundGraph)
        printSolution(g_before, upperBoundGraph)
        println("#MinCost == upperBound")
        println("#Finished: so")
        return upperBound
    else
        for op in reverse(path)
            if op.type == deleteOp
                setWeight(g, op.u, op.v, -Inf)
            elseif op.type == mergeOp
                merge!(g, op.u, op.v)
            elseif op.type == setForbiddenOp
                setForbidden(g, op.u, op.v)
            end
        end
    end

    g.foundSolution = true
    undoAll!(g)

    printSolution(g_before, g)
    println("#Finished: so")
    return minCost
end

function isConnected(g::Graph)::Bool
    visited = fill(false, g.n_total)
    q = Queue{Int}()

    cur = g.indexes[1]
    @inbounds visited[cur] = true
    enqueue!(q, cur)

    while !isempty(q)
        u = dequeue!(q)

        for v in g.indexes
            if hasEdge(g, u, v)
                if !visited[v]
                    @inbounds visited[v] = true
                    enqueue!(q, v)
                end
            end
        end
    end

    for cur = g.indexes
        @inbounds if !visited[cur] return false end
    end

    return true
end

function getComponentGraphs(g::Graph)::Array{Graph, 1}
    visited = fill(false, g.n_total)
    q = Queue{Int}()

    graphs::Array{Graph, 1} = []
    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        g_cur = Graph(g.n_total, g.weights, g.changeStack, [], fill(false, g.n_total), false)

        while !isempty(q)
            u = dequeue!(q)

            g_cur.indexes_lookup[u] = true
            push!(g_cur.indexes, u)

            for v in g.indexes
                if hasEdge(g, u, v)
                    if !visited[v]
                        @inbounds visited[v] = true
                        enqueue!(q, v)
                    end
                end
            end
        end

        push!(graphs, g_cur)
    end

    return graphs
end

function solveAll(g::Graph, params::Params)
    removeClusters(g)
    graphs = getComponentGraphs(g)

    total_k = 0
    for graph in graphs
        total_k += solveInstance(graph, params)
    end

    println("#last-k: $total_k")
    println("#recursive steps: $(recursive_steps)")
end


function parseParams()
    applyLB2 = false
    applyLB2Freq = 1

    applyLB3 = true
    applyLB3Freq = 2

    applyComponentSplitting = true
    applyComponentSplittingFreq = 7

    applyBranchDeleteFirst = true

    heuristicOuterSteps = 40
    heuristicInnerSteps = 50000

    applyMergeGreaterK = true
    applyMergeGreaterKFreq = 50

    applyHeavyNonEdge = true
    applyHeavyNonEdgeFreq = 40

    applyHeavyBothEnds = false
    applyHeavyBothEndsFreq = 3

    applyHeavySingleEnd = true
    applyHeavySingleEndFreq = 30

    applyLN = false
    applyLNFreq = 10

    for i in 5:size(ARGS, 1)
        if ARGS[i] == "-applyLB2Freq"
            applyLB2Freq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyLB3Freq"
            applyLB3Freq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyComponentSplittingFreq"
            applyComponentSplittingFreq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyBranchDeleteFirst"
            applyBranchDeleteFirst = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-heuristicOuterSteps"
            heuristicOuterSteps = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-heuristicInnerSteps"
            heuristicInnerSteps = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyLB2"
            applyLB2 = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyLB3"
            applyLB3 = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyComponentSplitting"
            applyComponentSplitting = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyMergeGreaterK"
            applyMergeGreaterK = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyMergeGreaterKFreq"
            applyMergeGreaterKFreq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyHeavyNonEdge"
            applyHeavyNonEdge = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyHeavyNonEdgeFreq"
            applyHeavyNonEdgeFreq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyHeavyBothEnds"
            applyHeavyBothEnds = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyHeavyBothEndsFreq"
            applyHeavyBothEndsFreq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyHeavySingleEnd"
            applyHeavySingleEnd = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyHeavySingleEndFreq"
            applyHeavySingleEndFreq = parse(Int, ARGS[i+1])
        elseif ARGS[i] == "-applyLN"
            applyLN = parse(Bool, ARGS[i+1])
        elseif ARGS[i] == "-applyLNFreq"
            applyLNFreq = parse(Int, ARGS[i+1])
        end
    end

    return Params(
        applyLB2,
        applyLB2Freq,
        applyLB3,
        applyLB3Freq,
        applyComponentSplitting,
        applyComponentSplittingFreq,
        applyBranchDeleteFirst,
        heuristicOuterSteps,
        heuristicInnerSteps,
        applyMergeGreaterK,
        applyMergeGreaterKFreq,
        applyHeavyNonEdge,
        applyHeavyNonEdgeFreq,
        applyHeavyBothEnds,
        applyHeavyBothEndsFreq,
        applyHeavySingleEnd,
        applyHeavySingleEndFreq,
        applyLN,
        applyLNFreq
    )
end

function main()
    params = parseParams()

    g = readMatrixFromInput()
    setupLogger()
    solveAll(g, params)
end

main()
