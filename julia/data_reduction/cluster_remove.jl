function removeClusters(g::Graph)
    cluster_components = []
    visited = fill(false, g.n_total)
    q = Queue{Int}()
    n = size(g.indexes, 1)

    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        isClusterComponent = true
        prev_degree = -1
        component_size = 0
        while !isempty(q)
            u = dequeue!(q)
            component_size += 1

            degree = 0
            for v = g.indexes
                if hasEdge(g, u, v)
                    degree += 1
                    if !visited[v]
                        @inbounds visited[v] = true
                        enqueue!(q, v)
                    end
                end
            end

            if prev_degree != -1 && prev_degree != degree
                isClusterComponent = false
            end

            prev_degree = degree
        end

        if component_size != 1 && component_size - 1 != prev_degree
            isClusterComponent = false
        end

        if isClusterComponent
            push!(cluster_components, cur)
        end
    end

    if isempty(cluster_components) return 0 end

    nodes_to_remove::IntArray = []
    for cur in cluster_components
        push!(nodes_to_remove, cur)
        for i = g.indexes
            if (hasEdge(g, cur, i))
                push!(nodes_to_remove, i)
            end
        end
    end

    # TODO: Maybe use set difference
    filter!(x -> !(x in nodes_to_remove), g.indexes)

    push!(g.changeStack, ClusterRemove(nodes_to_remove))
    return 1
end