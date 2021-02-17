using DataStructures
include("../helper/graph.jl")


function minimumCutPhase(g::Graph, N::Array{Int,1}, startNode::Int)
    A = [startNode]

    h = MutableBinaryMaxHeap{Float64}()

    handles = fill(-1, size(N, 1))
    for u = N
        if u == startNode continue end

        handle = push!(h, max(0, getWeight(g, startNode, u)))

        handles[handle] = u
    end

    while !isempty(h)
        value, handle = top_with_handle(h)
        addedNode = handles[handle]

        # Remove cur
        delete!(h, handle)

        # UpdateRemaining
        for node = h.nodes
            otherNode = handles[node.handle]
            if hasEdge(g, otherNode, addedNode)
                val = node.value + getWeight(g, otherNode, addedNode)
                update!(h, node.handle, val)
            end
        end

        push!(A, addedNode)
    end

    t = pop!(A)

    w = 0
    for u = A
        if g.indexes_lookup[u] && hasEdge(g, u, t)
            w += getWeight(g, u, t)
        end
    end


    s = last(A)
    mergeNonForbidden(g, s, t)

    return w, A
end

function minCutWeight(g::Graph, N::Array{Int,1})

    if size(N, 1) <= 1
        return 0
    end
    # cpy = copy(g)
    startNode = N[1]
    min_w = Inf
    merges = size(N, 1) - 1
    for i = 1:merges
        w, N = minimumCutPhase(g, N, startNode)
        min_w = min(w, min_w)
    end
    undoK(g, merges)

    # @assert graphEqual(g, cpy) "Not equal"


    return min_w
end

function test()
    g = readMatrixFromFile("/home/hjal/Desktop/Uni/AlgEng/julia/helper/tests/test_graphs/min_cut.dimacs")

    printMatrix(g)

    print(minCutWeight(g, [1, 2, 5, 6]))
    print(g.indexes)
end
