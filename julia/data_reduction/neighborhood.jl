include("../helper/graph.jl")
include("min_cut.jl")

function large_neighborhood(g::Graph, _::Float64, _::FloatArray, _::FloatArray)::Tuple{Int,Float64}
    for u in g.indexes
        N::IntArray = [u]
        rest::IntArray = []
        def::Float64 = 0
        cut::Float64 = 0

        for i in g.indexes
            if i == u continue end

            if hasEdge(g, u, i)
                push!(N, i)
            else
                push!(rest, i)
            end
        end
        n_neigh = size(N, 1)

        if n_neigh == 1 continue end

        found_zero_edge = false
        for i = 1:n_neigh
            @inbounds w = N[i]
            for j = i + 1:n_neigh
                @inbounds v = N[j]

                if g.weights[v, w] == 0.
                    found_zero_edge = true
                    break
                end

                if !hasEdge(g, v, w)
                    def += abs(getWeight(g, v, w))
                end
            end

            if found_zero_edge
                break
            end
        end

        if found_zero_edge
            continue
        end

        for i in N, j in rest
            if hasEdge(g, i, j)
                cut += getWeight(g, i, j)
            end
        end

        if 2 * def + cut < size(N, 1)
            k_decr = 0
            for v in N
                if v == u continue end

                if !hasEdge(g, u, v)
                    push!(path, Operation(addOp, u, v))
                    w = getWeight(g, u, v)
                    k_decr += abs(w)

                    push!(g.changeStack, AddEdge(u, v, w))

                    setWeight(g, u, v, Inf)
                end

                k_decr += merge!(g, u, v)
            end

            return size(N, 1) - 1, k_decr
        end

        minCut = minCutWeight(g, N)

        # println("mincut $(def + cut) $minCut $N")
        if def + cut < minCut
            k_decr = 0
            for v in N
                if v == u continue end

                if !hasEdge(g, u, v)
                    push!(path, Operation(addOp, u, v))
                    w = getWeight(g, u, v)
                    k_decr += abs(w)

                    push!(g.changeStack, AddEdge(u, v, w))

                    setWeight(g, u, v, Inf)
                end

                k_decr += merge!(g, u, v)
            end

            # println("LN $(size(N, 1) - 1), $k_decr, $N")
            return size(N, 1) - 1, k_decr
        end
    end

    return 0, 0
end

function large_neighborhood_rec(g::Graph)
    path::Array{Operation, 1} = []
    k_decr = 0.

    while true
        found = false
        for u in g.indexes
            N::IntArray = [u]
            rest::IntArray = []
            def::Float64 = 0
            cut::Float64 = 0

            for i in g.indexes
                if i == u continue end

                if hasEdge(g, u, i)
                    push!(N, i)
                else
                    push!(rest, i)
                end
            end
            n_neigh = size(N, 1)

            if n_neigh == 1 continue end

            found_zero_edge = false
            for i = 1:n_neigh
                @inbounds w = N[i]
                for j = i + 1:n_neigh
                    @inbounds v = N[j]

                    if g.weights[v, w] == 0.
                        found_zero_edge = true
                        break
                    end

                    if !hasEdge(g, v, w)
                        def += abs(getWeight(g, v, w))
                    end
                end

                if found_zero_edge
                    break
                end
            end

            if found_zero_edge
                continue
            end

            for i in N, j in rest
                if hasEdge(g, i, j)
                    cut += getWeight(g, i, j)
                end
            end

            if 2 * def + cut < size(N, 1)
                for v in N
                    if v == u continue end

                    if !hasEdge(g, u, v)
                        push!(path, Operation(addOp, u, v))
                        w = getWeight(g, u, v)
                        k_decr += abs(w)

                        push!(g.changeStack, AddEdge(u, v, w))

                        setWeight(g, u, v, Inf)
                    end

                    k_decr += merge!(g, u, v)
                    push!(path, Operation(mergeOp, u, v))
                end
                found = true
            end

            if found break end

            minCut = minCutWeight(g, N)

            if def + cut < minCut
                for v in N
                    if v == u continue end

                    if !hasEdge(g, u, v)
                        push!(path, Operation(addOp, u, v))
                        w = getWeight(g, u, v)
                        k_decr += abs(w)

                        push!(g.changeStack, AddEdge(u, v, w))

                        setWeight(g, u, v, Inf)
                    end

                    k_decr += merge!(g, u, v)
                    push!(path, Operation(mergeOp, u, v))
                end
                found = true
            end

            if found break end
        end

        if !found
            break
        end
    end

    return path, k_decr
end