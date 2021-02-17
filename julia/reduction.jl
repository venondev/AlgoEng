include("./helper/graph.jl")
include("./helper/p3.jl")
include("./data_reduction/data_reduction.jl")

function addEdge(g, u, v, prev)
    return abs(getWeight(g, u, v)) + prev
end

function main()
    g = readMatrixFromInput()
    max = loopEdges(g, addEdge, 0.)

    k_decr = doDataReductionPre(g)

    # Remove cluster
    findP3AdjList(g)

    sort!(g.indexes)

    println("$(size(g.indexes, 1))")
    indexes = []
    i = 1
    for j = g.indexes
        push!(indexes, (j, i))
        i += 1
    end

    for i = 1:size(indexes, 1)
        u, u_i = indexes[i]
        for j = i+1:size(indexes, 1)
            v, v_i = indexes[j]
            w = getWeight(g, u, v)

            if abs(w) == Inf
                println("$u_i $v_i $(Int(sign(w) * max))")
            else
                println("$u_i $v_i $(Int(w))")
            end
        end
    end
    println("#weight: $k_decr")
end

main()
