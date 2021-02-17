include("../main.jl")

path = "/home/hjal/CLionProjects/AlgoEng/students"

function checkFor(p)
    g = readMatrixFromFile(p)
    k = solve(g)
    open("$p.solution", "r") do io
        solFile = read(io, String)
        n = parse(Int, solFile)

        if n != k
            println("$p $n $k")
            exit()
        else
            println("$p \xE2\x9C\x94")
        end
    end
end

checkFor("/home/hjal/CLionProjects/AlgoEng/students/1-random/r007.dimacs")

function run()
    for (root, dirs, files) in walkdir(path)
        println("Checking $root")
        for f = files
            if endswith(f, ".dimacs")
                checkFor("$root/$f")
            end
        end
    end
end

run()