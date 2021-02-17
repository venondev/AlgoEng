include("../helper/graph.jl")

path = ARGS[1]

function loopEdgesAcc(g::Graph, func, state)
    n = size(g.indexes, 1)
    for i = 1:n
        w = g.indexes[i]
        for j = i+1:n
            v = g.indexes[j]
            state = func(g, v, w, state)
        end
    end
    return state
end

function countRealNumEdges(g, u, v, prev)
    if abs(getWeight(g, u, v)) != Inf && getWeight(g, u, v) > 0
        return prev + 1
    end
    return prev
end

if isfile(path)
    g = readMatrixFromFile(path)
    println(loopEdgesAcc(g, countRealNumEdges, 0))
end
