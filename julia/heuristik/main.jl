include("../helper/graph.jl")
include("../helper/p3.jl")
include("../data_reduction/data_reduction.jl")

Base.exit_on_sigint(false)

struct ClusterGraph
    graph::Graph
    comps::Array{Array{Int, 1}, 1}
end

function getRandomCompIdx(clusterGraph::ClusterGraph, forbidden::Int)::Int
    comps = size(clusterGraph.comps, 1)
    if forbidden != -1
        if forbidden == 1
            return rand(forbidden+1:comps)
        elseif forbidden == comps
            return rand(1:comps-1)
        else
            return rand() > (forbidden / comps) ? rand(forbidden + 1:comps) : rand(1:forbidden - 1) 
        end
    else
        return rand(1:comps)
    end
end

P3List = Array{P3,1}

function lowerBound1(g::Graph, existingP3s::P3List, start::Int)::Int
    counter = fill(0, (g.n_total, g.n_total))

    total_p3s = 0
    s = size(existingP3s, 1)
    for edge_idx = start:s
        @inbounds p = existingP3s[edge_idx]
        if @inbounds g.indexes_lookup[p[1]] &&
            @inbounds g.indexes_lookup[p[2]] &&
            @inbounds g.indexes_lookup[p[3]] &&
            checkP3(g, p)
            total_p3s += 1

            s, m, e = p

            counter[s, m] += 1
            counter[m, s] += 1

            counter[m, e] += 1
            counter[e, m] += 1

            counter[s, e] += 1
            counter[e, s] += 1
        end
    end

    minValue = Inf
    u = -1
    v = -1
    for i in g.indexes, j in g.indexes
        if counter[i, j] != 0 && getWeight(g, i, j) != -Inf
            temp = abs(getWeight(g, i, j)) / counter[i, j]
            if temp < minValue
                minValue = temp
                u = i
                v = j
            end
        end
    end

    #println("minValue: $minValue")

    if hasEdge(g, u, v)
        setWeight(g, u, v, -Inf)
    else
        setWeight(g, u, v, Inf)
        merge!(g, u, v)
    end


    return total_p3s - counter[u, v]
end


function getRandomClusterGraph(g::Graph)::ClusterGraph
    clusterGraph = ClusterGraph(createGraph(fill(-1, (g.n_total, g.n_total))), [])
    nodes_to_process = copy(g.indexes)

    while !isempty(nodes_to_process)
        cur = nodes_to_process[rand(1:size(nodes_to_process, 1))]

        cluster = [cur]
        for i = g.indexes
            if hasEdge(g, cur, i) && i in nodes_to_process
                push!(cluster, i)
            end
        end

        for u = cluster, v = cluster
            setWeight(clusterGraph.graph, u, v, Inf)
        end

        push!(clusterGraph.comps, cluster)

        filter!(x -> !(x in cluster), nodes_to_process)
    end

    return clusterGraph
end

function heavyNonEdge!(g::Graph, N::FloatArray, g_cluster::Graph)::Float64
    total_saved = 0.
    for i = g.indexes, j = g.indexes
        if !hasEdge(g, i, j) && hasEdge(g_cluster, i, j) && abs(getWeight(g, i, j)) > @inbounds N[j]
            t = 0
            for v = g.indexes
                if hasEdge(g, v, j)
                    t += getWeight(g, v, j)
                end
                setWeight(g_cluster, v, j, -Inf)
            end
            w = abs(getWeight(g, i, j))
            total_saved += w - t
        end
    end
    return total_saved
end

function splitComp(cost::Float64, before::Graph, clusterGraph::ClusterGraph)::Float64
    c_idx = getRandomCompIdx(clusterGraph, -1)
    c = clusterGraph.comps[c_idx]

    if size(c, 1) == 1 return cost end

    after_cost = cost

    splitat = rand(1:size(c, 1) - 1)
    for u = c[1:splitat], v = c[splitat+1:end]
        if hasEdge(before, u, v)
            after_cost += abs(getWeight(before, u, v))
        else
            after_cost -= abs(getWeight(before, u, v))
        end
    end

    if after_cost < cost
        #println("Unmerge $after_cost")
        for u = c[1:splitat], v = c[splitat+1:end]
            setWeight(clusterGraph.graph, u, v, -Inf)
        end

        clusterGraph.comps[c_idx] = c[1:splitat]
        push!(clusterGraph.comps, c[splitat+1:end])
    end

    return min(cost, after_cost)
end

function mergeComps(cost::Float64, before::Graph, clusterGraph::ClusterGraph)::Float64
    if size(clusterGraph.comps, 1) == 1 return cost end

    c1_idx = getRandomCompIdx(clusterGraph, -1)
    c2_idx = getRandomCompIdx(clusterGraph, c1_idx)

    c1 = clusterGraph.comps[c1_idx]
    c2 = clusterGraph.comps[c2_idx]

    after_cost = cost

    for u in c1, v in c2
        if hasEdge(before, u, v)
            after_cost -= abs(getWeight(before, u, v))
        else
            after_cost += abs(getWeight(before, u, v))
        end
    end

    if after_cost < cost
        #println("Merge $after_cost")
        for u in c1, v in c2
            setWeight(clusterGraph.graph, u, v, Inf)
        end

        clusterGraph.comps[c1_idx] = vcat(c1, c2)
        deleteat!(clusterGraph.comps, c2_idx)
    end

    return min(cost, after_cost)
end

function moveVertexToOtherComp(cost::Float64, before::Graph, clusterGraph::ClusterGraph)
    if size(clusterGraph.comps, 1) == 1 return cost end

    c1_idx = getRandomCompIdx(clusterGraph, -1)
    c2_idx = getRandomCompIdx(clusterGraph, c1_idx)

    c1 = clusterGraph.comps[c1_idx]
    c2 = clusterGraph.comps[c2_idx]

    v_idx = rand(1:size(c1, 1))
    vToMove = c1[v_idx]

    after_cost = cost

    for v in c1
        if v == vToMove continue end

        if hasEdge(before, vToMove, v)
            after_cost += abs(getWeight(before, vToMove, v))
        else
            after_cost -= abs(getWeight(before, vToMove, v))
        end
    end

    for v in c2
        if hasEdge(before, vToMove, v)
            after_cost -= abs(getWeight(before, vToMove, v))
        else
            after_cost += abs(getWeight(before, vToMove, v))
        end
    end

    if after_cost < cost
        #println("Move $after_cost")
        for v in c1
            if v == vToMove continue end

            setWeight(clusterGraph.graph, vToMove, v, -Inf)
        end

        for v in c2
            setWeight(clusterGraph.graph, vToMove, v, Inf)
        end

        push!(c2, vToMove)
        deleteat!(c1, v_idx)

        if isempty(c1)
            deleteat!(clusterGraph.comps, c1_idx)
        end
    end

    return min(cost, after_cost)
end

function solveRandom(before::Graph, maxSteps::Int)::Tuple{ClusterGraph, Float64}
    clusterGraph = getRandomClusterGraph(before)

    cost = getSolutionSize(before, clusterGraph.graph)

    for i = 1:maxSteps
        step = rand(1:3)

        if step == 1
            cost = mergeComps(cost, before, clusterGraph)
        elseif step == 2
            cost = moveVertexToOtherComp(cost, before, clusterGraph)
        else 
            cost = splitComp(cost, before, clusterGraph)
        end
    end

    return clusterGraph, cost
end

function solveHeuristik(before::Graph, params::Params)::Tuple{Graph, Float64}
    minCluster = ClusterGraph(createGraph(fill(-1, (1, 1))), [])
    minCost = Inf
    for i = 1:params.heuristicOuterSteps
        clusterGraphTemp, costTemp = solveRandom(before, params.heuristicInnerSteps)
        if costTemp < minCost
            minCluster = clusterGraphTemp
            minCost = costTemp
        end
    end
    return minCluster.graph, minCost
end

function main()
    before = readMatrixFromInput()
    solveHeuristik(before)
end