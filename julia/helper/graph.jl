using DataStructures
include("./math.jl")

# Types
P3 = Tuple{Int,Int,Int}
FloatArray = Array{Float64,1}
FloatMatrix = Array{Float64,2}
IntArray = Array{Int,1}
BoolArray = Array{Bool,1}

struct Merge
    u_weights::FloatArray
    u::Int
    v::Int
    k_decr::Float64
end

struct ClusterRemove
    indexes::Array{Int64,1}
end

struct AddEdge
    u::Int
    v::Int
    weight::Float64
end

struct SetForbidden
    u::Int
    v::Int
    weight::Float64
end
struct MergeNonForbidden
    u_weights::FloatArray
    u::Int
    v::Int
end

Change = Union{Merge,MergeNonForbidden,ClusterRemove,SetForbidden, AddEdge}
mutable struct Graph
    n_total::Int
    weights::FloatMatrix
    changeStack::Stack{Change}
    indexes::Array{Int64,1}
    indexes_lookup::Array{Bool,1}
    foundSolution::Bool
end

@enum OpType begin
    mergeOp = 1
    deleteOp = 2
    setForbiddenOp = 3
    addOp = 4
end

struct Operation
    type::OpType
    u::Int
    v::Int
end

struct Params
    applyLB2::Bool
    applyLB2Freq::Int
    applyLB3::Bool
    applyLB3Freq::Int
    applyComponentSplitting::Bool
    applyComponentSplittingFreq::Int
    applyBranchDeleteFirst::Bool
    heuristicOuterSteps::Int
    heuristicInnerSteps::Int
    applyMergeGreaterK::Bool
    applyMergeGreaterKFreq::Int
    applyHeavyNonEdge::Bool
    applyHeavyNonEdgeFreq::Int
    applyHeavyBothEnds::Bool
    applyHeavyBothEndsFreq::Int
    applyHeavySingleEnd::Bool
    applyHeavySingleEndFreq::Int
    applyLN::Bool
    applyLNFreq::Int
end

function createGraph(matrix)
    n = size(matrix, 1)

    g = Graph(n, matrix, Stack{Change}(), [i for i in 1:n], fill(true, n), false)
    return g
end

function graphEqual(g::Graph, g_comp::Graph)
    sort!(g.indexes)
    sort!(g_comp.indexes)

    for i = 1:size(g.indexes, 1)
        if g.indexes[i] != g_comp.indexes[i]
            println("Compare Index error, $(g.indexes) $(g_comp.indexes)")
            return false
        end
    end

    for i in g.indexes, j in g.indexes
        if g.weights[i, j] != g_comp.weights[i, j]
            println("Weights error $i $j: $(g.weights[i, j]) != $(g_comp.weights[i, j])")
            return false
        end
    end

    return true
end

function getWeight(g::Graph, u::Int, v::Int)::Float64
    return @inbounds g.weights[u, v]
end

function setWeight(g::Graph, u::Int, v::Int, value::Float64)
    @inbounds g.weights[u, v] = value
    @inbounds g.weights[v, u] = value
end

function hasEdge(g::Graph, u::Int, v::Int)::Bool
    return @inbounds g.weights[u, v] > 0
end

function merge!(g::Graph, u::Int, v::Int)::Float64

    u_weights = fill(-1.0, g.n_total)
    k_decr = 0.
    filter!(x -> x != v, g.indexes)
    @inbounds g.indexes_lookup[v] = false

    for i = g.indexes
        if i == u || i == v continue end

        u_i = getWeight(g, u, i)
        v_i = getWeight(g, v, i)

        @inbounds u_weights[i] = u_i

        # decrease k and update has_edge
        if hasEdge(g, u, i) != hasEdge(g, v, i)
            k_decr += min(abs(u_i), abs(v_i))
        end

        # update weights
        setWeight(g, u, i, u_i + v_i)
    end

    push!(g.changeStack, Merge(u_weights, u, v, k_decr))
    return k_decr
end

function mergeNonForbidden(g::Graph, u::Int, v::Int)
    u_weights = fill(-1.0, g.n_total)
    filter!(x -> x != v, g.indexes)
    @inbounds g.indexes_lookup[v] = false

    for i = g.indexes
        if i == u || i == v continue end

        u_weights[i] = getWeight(g, u, i)

        u_i = max(0, getWeight(g, u, i))
        v_i = max(0, getWeight(g, v, i))

        setWeight(g, u, i, u_i + v_i)
    end

    push!(g.changeStack, MergeNonForbidden(u_weights, u, v))
end

function setForbidden(g::Graph, u::Int, v::Int)
    forbidden = SetForbidden(u, v, getWeight(g, u, v))
    setWeight(g, u, v, forbidden.weight > 0 ? Inf : -Inf)
    push!(g.changeStack, forbidden)
end

function undoAll!(g::Graph)
    while !isempty(g.changeStack)
        undo(g, pop!(g.changeStack))
    end
end

function undoK(g::Graph, k::Int)
    for i = 1:k
        if isempty(g.changeStack) break end
        undo(g, pop!(g.changeStack))
    end
end

function undo(g::Graph, addEdge::AddEdge)
    setWeight(g, addEdge.u, addEdge.v, addEdge.weight)
end

function undo(g::Graph, merge::Merge)
    for i = g.indexes
        if i == merge.u || i == merge.v continue end

        if g.foundSolution
            setWeight(g, merge.v, i, getWeight(g, merge.u, i))
        else
            setWeight(g, merge.u, i, @inbounds merge.u_weights[i])
        end
    end

    push!(g.indexes, merge.v)
    @inbounds  g.indexes_lookup[merge.v] = true
end

function undo(g::Graph, merge::MergeNonForbidden)
    for i = g.indexes
        if i == merge.u || i == merge.v continue end

        setWeight(g, merge.u, i, @inbounds merge.u_weights[i])
    end

    push!(g.indexes, merge.v)
    @inbounds  g.indexes_lookup[merge.v] = true
end

function undo(g::Graph, clusterRemove::ClusterRemove)
    append!(g.indexes, clusterRemove.indexes)
end

function undo(g::Graph, forbidden::SetForbidden)
    setWeight(g, forbidden.u, forbidden.v, forbidden.weight)
end

function getSolutionSize(before::Graph, after::Graph)::Float64
    w = 0

    for i = 1:size(before.indexes, 1)
        u = before.indexes[i]
        for j = i + 1:size(before.indexes, 1)
            v = before.indexes[j]
            if (hasEdge(before, u, v) && !hasEdge(after, u, v)) ||
               (!hasEdge(before, u, v) && hasEdge(after, u, v))
                w += abs(getWeight(before, u, v))
            end
        end
    end

    return w
end

function printRows(g::Graph, rows::IntArray)
    print("\t")
       
    for j = 1:g.n_total
        print("  $j:\t")
    end
    println(" ")

    for i = rows
        print("$i:\t")
        for j = 1:g.n_total
            if !g.indexes_lookup[j]
                print(" _\t")
            else
                print(" $(g.weights[i, j])\t")
            end
        end
        println("")
    end
end

function printMatrix(g::Graph, showNonEdge::Bool=false)
    print("\t")
    for i = g.indexes
        print("  $i:\t")
    end
    println(" ")

    for i = g.indexes
        print("$i:\t")
        for j = g.indexes
            if showNonEdge
                print(" $(g.weights[i, j])\t")
            else 
                if hasEdge(g, i, j)
                    print(" $(g.weights[i, j])\t")
                else
                    print(" _\t")
                end
            end
        end
        println("")
    end
end

function printEdges(g::Graph)
    println("========================")
    n = size(g.indexes, 1)

    for i = g.indexes
        println("$i")
    end

    #println("$n")


    for i = 1:size(g.indexes, 1)
        u = g.indexes[i]
        for j = i + 1:size(g.indexes, 1)
            v = g.indexes[j]
            if g.weights[u, v] > 0
                println("$u $v $(g.weights[u, v])")
            end
        end
    end
    println("========================")
end

function printAdjList(g::Graph)
    list = matrixToList(g)
    for i = 1:g.n_total
        println("$i : $(list[i])")
    end
end

function printSolution(before::Graph, after::Graph)
    for i = 1:size(before.indexes, 1)
        u = before.indexes[i]
        for j = i + 1:size(before.indexes, 1)
            v = before.indexes[j]
            if (hasEdge(before, u, v) && !hasEdge(after, u, v)) ||
               (!hasEdge(before, u, v) && hasEdge(after, u, v))
                println("$u $v")
            end
        end
    end
end

function loopEdges(g::Graph, func, start)
    n = size(g.indexes, 1)
    for i = 1:n
        w = g.indexes[i]
        for j = i + 1:n
            v = g.indexes[j]
            start = func(g, v, w, start)
        end
    end
    return start
end

function readMatrixFromFile(path)::Graph
    open(path, "r") do io
        f = read(io, String)

        lines = split(f, "\n")
        pop!(lines)
        n = parse(Int, popat!(lines, 1))

        matrix = fill(-1., (n, n))

        for line = lines
            pieces = split(line, ' ', keepempty=false)
            pieces = map(x -> parse(Int, x), pieces)

            s, e, weight = pieces
            matrix[s, e] = weight
            matrix[e, s] = weight
        end
        return createGraph(matrix)
    end
end

function readMatrixFromInput()::Graph
    n = parse(Int, readline())

    matrix = fill(-1., (n, n))

    lines = bino(n, 2)
    for i = 0:lines
        line = readline();

        if line == ""
            break
        end

        pieces = split(line, ' ', keepempty=false)
        pieces = map(x -> parse(Int, x), pieces)

        s, e, weight = pieces
        matrix[s, e] = weight
        matrix[e, s] = weight
    end

    return createGraph(matrix)
end

function readSolutionFile(path)
    try
        open("$path.solution", "r") do io
            solution = parse(Int, read(io, String))
            return solution
        end
    catch
        return -1
    end
end

function copy(m::Merge)
    return Merge(copy(m.u_weights), m.u, m.v, m.k_decr)
end

function copy(s::SetForbidden)
    return SetForbidden(s.u, s.v, s.weight)
end

function copy(c::ClusterRemove)
    return ClusterRemove(copy(c.indexes))
end

function copy(m::FloatArray)
    ret = fill(-1., size(m, 1))

    for i = 1:size(ret, 1)
        ret[i] = m[i]
    end

    return ret
end

function copy(m::IntArray)
    ret = fill(5, size(m, 1))

    for i = 1:size(ret, 1)
        ret[i] = m[i]
    end

    return ret
end

function copy(m::BoolArray)
    ret = fill(false, size(m, 1))

    for i = 1:size(ret, 1)
        ret[i] = m[i]
    end

    return ret
end

function copy(m::FloatMatrix)
    ret = fill(-1.0, size(m))

    for i in 1:size(ret, 1), j in 1:size(ret, 2)
        ret[i, j] = m[i, j]
    end

    return ret
end

function copy(s::Stack{Change})
    ret = Stack{Change}()
    ret2 = Stack{Change}()
    for m in s
        push!(ret, copy(m))
    end

    for m in ret
        push!(ret2, m)
    end

    return ret2
end

function copy(g::Graph)
    return Graph(g.n_total, copy(g.weights), copy(g.changeStack), copy(g.indexes), copy(g.indexes_lookup), g.foundSolution)
end
