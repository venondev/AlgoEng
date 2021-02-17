include("../helper/graph.jl")
include("../data_reduction/data_reduction.jl")
include("../helper/p3.jl")

path = "/home/hjal/CLionProjects/AlgoEng/students"
save_path = "/home/hjal/CLionProjects/AlgoEng/benchmark/data_reduction/"
println(path)

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
    if abs(getWeight(g, u, v)) != Inf
        return prev + 1
    end
    return prev
end

file_name = ARGS[1]
description = ARGS[2]

open("$save_path/$file_name.txt", "w") do io
    write(io, description)
end

open("$save_path/$file_name.csv", "w") do io
    header = "file,n_pre,n_post,e_pre,e_post,d,time\n"

    write(io, header)
    println(header)

    for (root, dirs, files) in walkdir(path)
        println("Checking $root")
        for f = files
            if endswith(f, ".dimacs")
                g = readMatrixFromFile("$root/$f")

                e_pre = loopEdgesAcc(g, countRealNumEdges, 0)
                n_pre = size(g.indexes, 1)

                d = 0
                time = @elapsed begin
                    # Data reduction rules
                    d = doDataReductionPre(g)

                    # # Remove cluster
                    findP3AdjList(g)
                end

                e_post = loopEdgesAcc(g, countRealNumEdges, 0)
                n_post = size(g.indexes, 1)

                line = "$f,$n_pre,$n_post,$e_pre,$e_post,$d,$time\n"

                write(io, line)
                print(line)
            end
        end
    end
end
