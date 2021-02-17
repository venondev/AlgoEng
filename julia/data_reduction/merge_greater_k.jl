include("../helper/graph.jl")

function mergeGreaterK(g::Graph, k::Float64)::Tuple{Int, Float64}
    if k > 0  
        for i = 1:size(g.indexes, 1)
            @inbounds u = g.indexes[i]
            for j = i:size(g.indexes, 1)
                @inbounds v = g.indexes[j]
                if hasEdge(g, v, u) && getWeight(g, v, u) > k
                    return 1, merge!(g, u, v)
                end
            end
        end
    end

    return nothing
end

function applyMergeGreaterKRec(g::Graph, k::Float64)
    path::Array{Operation, 1} = []
    k_decr = 0.
    while true
        found = false
        if (k - k_decr) > 0  
            for i = 1:size(g.indexes, 1)
                @inbounds u = g.indexes[i]
                for j = i:size(g.indexes, 1)
                    @inbounds v = g.indexes[j]
                    if hasEdge(g, v, u) && getWeight(g, v, u) > k
                        k_decr += merge!(g, u, v)
                        push!(path, Operation(mergeOp, u, v))
                        found = true
                        break
                    end
                end
                if found
                    break
                end
            end
        end

        if !found
            break
        end
    end

    return path, k_decr
end
