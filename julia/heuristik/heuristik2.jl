include("../helper/graph.jl")
include("../helper/p3.jl")

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

    if hasEdge(g, u, v)
        setWeight(g, u, v, -Inf)
    else
        setWeight(g, u, v, Inf)
        merge!(g, u, v)
    end


    return total_p3s - counter[u, v]
end

function printSolution(before::Graph, after::Graph)
    w = 0
    for i = 1:before.n_total
        for j = i + 1:before.n_total
            if (hasEdge(before, i, j) && !hasEdge(after, i, j)) ||
               (!hasEdge(before, i, j) && hasEdge(after, i, j))
                println("$i $j")
                w += abs(getWeight(before, i, j))
            end
        end
    end

    println("#last-k: $w")
end

function main()
    g = readMatrixFromInput()
    before = copy(g)

    found, p3s = findP3AdjList(g)

    iter = 0
    while found
        initial_total_p3s = size(p3s, 1)
        while true
            total_p3s_left = lowerBound1(g, p3s, 1)

            if total_p3s_left / initial_total_p3s < 0.5
                break
            end
        end

        found, p3s = findP3AdjList(g)
        iter += 1
    end

    g.foundSolution = true
    undoAll!(g)

    printSolution(before, g)
end

main() 