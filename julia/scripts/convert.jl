include("../main.jl")

path = "/home/hjal/CLionProjects/AlgoEng/students/3-actionseq/a003.dimacs"

open(path, "r") do io
    println("Starting")
    g = readMatrixFromFile(path)
    println("read graph")
    solution = solve(g)

    f = read(io, String)

    lines = split(f, "\n")
    pop!(lines)
    n = parse(Int, popat!(lines, 1))

    println("==========================================")
    for i = 1:n println(i) end

    total_weight = 0
    for line in lines
        pieces = split(line, ' ', keepempty=false)
        pieces = map(x -> parse(Int, x), pieces)

        s, e, weight = pieces

        for (a, b) in solution
            if a == s && b == e
                weight *= -1
                total_weight += abs(weight)
                break
            end
        end


        if weight > 0
            println("$s $e")
        end
    end
    println("==========================================")
    println("#weight: $total_weight")
end
