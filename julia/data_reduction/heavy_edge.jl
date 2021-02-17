
function findAdjWeights(g::Graph)::Tuple{FloatArray, FloatArray}
    N::FloatArray = fill(0.0, g.n_total)
    V::FloatArray = fill(0.0, g.n_total)

    for u in g.indexes, v in g.indexes
        if u == v continue end
        w = getWeight(g, u, v)

        @inbounds V[u] += abs(w)

        if hasEdge(g, u, v)
            @inbounds N[u] += w
        end
    end

    return N, V
end

function heavyNonEdge!(g::Graph, N::FloatArray)::Int
    changes = 0
    for i = g.indexes, j = g.indexes
        if !hasEdge(g, i, j) && abs(getWeight(g, i, j)) > @inbounds N[j]
            setForbidden(g, i, j)
            changes += 1
        end
    end
    return changes
end

function heavyEdgeSingleEnd!(g::Graph, k::Float64, N::FloatArray, V::FloatArray)::Tuple{Int, Float64}
    for u = g.indexes, v = g.indexes
        if v == u continue end

        w = getWeight(g, u, v)
        if hasEdge(g, u, v) && w >= @inbounds V[u] - abs(w)
            return 1, merge!(g, u, v)
        end
    end

    return 0, 0.
end

function heavyEdgeBothEnds!(g::Graph, k::Float64, N::FloatArray, V::FloatArray)::Tuple{Int, Float64}
    for u = g.indexes, v = g.indexes
        if v == u continue end

        w = getWeight(g, u, v)

        if hasEdge(g, u, v) && w >= (@inbounds N[u] - w) + (@inbounds N[v] - w)
            return 1, merge!(g, u, v)
        end
    end

    return 0, 0.
end

function heavyEdgeEndsRec(g::Graph, both::Bool, single::Bool, forbidden::Bool)
    path::Array{Operation, 1} = []
    k_decr = 0.

    while true
        N, V = findAdjWeights(g)

        found = false
        for u = g.indexes, v = g.indexes
            if v == u continue end

            w = getWeight(g, u, v)
            if w == -Inf continue end

            if forbidden && !hasEdge(g, u, v) && abs(w) > @inbounds N[v]
                setForbidden(g, u, v)
                push!(path, Operation(setForbiddenOp, u, v))
            end

            if both && hasEdge(g, u, v) && w >= (@inbounds N[u] - w) + (@inbounds N[v] - w)
                k_decr += merge!(g, u, v)
                push!(path, Operation(mergeOp, u, v))
                found = true
                break
            end

            if single && hasEdge(g, u, v) && w >= @inbounds V[u] - abs(w)
                k_decr += merge!(g, u, v)
                push!(path, Operation(mergeOp, u, v))
                found = true
                break
            end
        end
        if !found break end
    end

    return path, k_decr
end
