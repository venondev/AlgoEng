
include("../helper/graph.jl")
include("../data_reduction/data_reduction.jl")

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
    g_before = copy(g)

    k_decr = doDataReductionPre(g)
    n = size(g.indexes, 1)
    println("# $(g.n_total) $n")

    if n >= 3
        graph_string = "$(g.n_total)\n"
        sort!(g.indexes)
        for i = 1:n
            u = g.indexes[i]
            for j = i + 1:n
                v = g.indexes[j]
                w = getWeight(g, u, v)

                graph_string = graph_string * ("$u $v $(w)\n")
            end
        end

        out = read(pipeline(`echo $graph_string`, `python3 ./julia/ortools/main.py`), String)

        for line in split(out, "\n")
            if line == "" continue end
            pieces = split(line, ' ', keepempty=false)
            pieces = map(x -> parse(Int, x), pieces)

            s, e, has_edge = pieces
            if has_edge == 1
                setWeight(g, s, e, Inf)
            else
                setWeight(g, s, e, -Inf)
            end
        end

    end

    g.foundSolution = true
    undoAll!(g)

    printSolution(g_before, g)
end

main()