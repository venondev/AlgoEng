include("../helper/graph.jl")

function gurobiLB(g::Graph)
    n = size(g.indexes, 1)
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

    out = read(pipeline(`echo $graph_string`, `python3 ./julia/gurobi/main.py`, `grep "Obj:"`), String)
    _, objValue = split(out, " ")
    return parse(Float64, objValue)
end

function main()
    before = readMatrixFromInput()
    println(gurobiLB(before))
end