include("graph.jl")

using DataStructures

function checkP3(g::Graph, (s, m, e)::P3)::Bool
    return hasEdge(g, s, m) && hasEdge(g, m, e) && !hasEdge(g, s, e) && s != e
end

function minWeight(g::Graph, (s, m, e)::P3)::Float64
    return min(getWeight(g, s, m), getWeight(g, m, e), abs(getWeight(g, s, e)))
end

function matrixToList(g::Graph)::Array{IntArray,1}
    ret::Array{Array{Int,1},1} = fill([], g.n_total)

    for i = 1:g.n_total ret[i] = [] end

    for i in g.indexes, j in g.indexes
        if hasEdge(g, i, j)
            push!(ret[i], j)
        end
    end

    return ret
end

function sortFunc(x::Tuple{P3,Float64})::Int
    p3, minWeight = x
    return minWeight
end

function bfsDepth2(g::Graph, start::Int, p3s::Set{Tuple{P3,Float64}})
    visited = fill(false, g.n_total)
    q = Queue{Int}()

    enqueue!(q, start)

    while !isempty(q)
        u = dequeue!(q)

        for v in g.indexes
            if hasEdge(g, u, v)
                if !visited[v] && u == start
                    @inbounds visited[v] = true
                    enqueue!(q, v)
                end

                if !hasEdge(g, start, v) && start != v
                    p = (start, u, v)
                    if start > v
                        p = (v, u, start)
                    end
                    
                    push!(p3s, (p, minWeight(g, p)))
                end
            end
        end
    end
end

function findP3sWithWeightFast(g::Graph)::Set{Tuple{P3, Float64}}
    visited = fill(false, g.n_total)
    q = Queue{Int}()

    degree = fill(0, g.n_total)
    ccsize = fill(0, g.n_total)
    k = 1
    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        compSize = 0
        while !isempty(q)
            compSize += 1
            u = dequeue!(q)

            degreeU = 0
            for v in g.indexes
                if hasEdge(g, u, v)
                    if !visited[v]
                        @inbounds visited[v] = true
                        enqueue!(q, v)
                    end
                    degreeU += 1
                end
            end
            degree[u] = degreeU
        end
        ccsize[cur] = compSize
    end

    visited = fill(false, g.n_total)
    q = Queue{Int}()

    p3s::Set{Tuple{P3,Float64}} = Set()
    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        while !isempty(q)
            u = dequeue!(q)

            if degree[u] < ccsize[cur] - 1
                # BFS Tiefe 2
                bfsDepth2(g, u, p3s)
            end

            for v in g.indexes
                if hasEdge(g, u, v)
                    if !visited[v]
                        @inbounds visited[v] = true
                        enqueue!(q, v)
                    end
                end
            end
        end
    end

    return p3s
end

function findP3sWithWeight(g::Graph)
    list = matrixToList(g)
    visited = fill(false, g.n_total)
    q = Queue{Int}()

    p3s::Array{Tuple{P3,Float64},1} = []
    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        while !isempty(q)
            m = dequeue!(q) # Middle Node 
            n_neigh = size(list[m], 1)

            for i = 1:n_neigh 
                @inbounds s = list[m][i] # Start node
                if !visited[s]
                    @inbounds visited[s] = true
                    enqueue!(q, s)
                end

                for j = i + 1:n_neigh 
                    @inbounds e = list[m][j] # End Node
                    p = (s, m, e)

                    if checkP3(g, p)
                        push!(p3s, (p, minWeight(g, p)))
                    end
                end
            end
        end
    end

    return p3s
end

function findP3sWithWeightSorted(g::Graph)::Array{P3,1}
    p3s = findP3sWithWeightFast(g)
    sort!(p3s, by=sortFunc, rev=true)
    return p3s
end

function findP3AdjList(g::Graph)::Tuple{Bool,Array{P3,1}}
    p3s = findP3sWithWeight(g)

    if isempty(p3s) return false, [] end

    p3sArray = collect(p3s)

    sort!(p3sArray, by=sortFunc, rev=true)

    p3s_conv = map(first, p3sArray)

    return true, p3s_conv 
end

function isClusterGraph(g::Graph)::Bool
    visited = fill(false, g.n_total)
    q = Queue{Int}()

    degree = fill(0, g.n_total)
    ccsize = fill(0, g.n_total)
    k = 1
    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        compSize = 0
        while !isempty(q)
            compSize += 1
            u = dequeue!(q)

            degreeU = 0
            for v in g.indexes
                if hasEdge(g, u, v)
                    if !visited[v]
                        @inbounds visited[v] = true
                        enqueue!(q, v)
                    end
                    degreeU += 1
                end
            end
            degree[u] = degreeU
        end
        ccsize[cur] = compSize
    end

    visited = fill(false, g.n_total)
    q = Queue{Int}()

    for cur = g.indexes
        @inbounds if visited[cur] continue end

        @inbounds visited[cur] = true
        enqueue!(q, cur)

        while !isempty(q)
            u = dequeue!(q)

            if degree[u] < ccsize[cur] - 1
                return false
            end
        end
    end

    return true
end