include("../helper/graph.jl")

path = "/home/lm/Desktop/Uni/AlgEng/students/3-actionseq"

for (root, dirs, files) in walkdir(path)
    println("Checking $root")
    for f = files
        if endswith(f, ".dimacs")
            g = readMatrixFromFile("$root/$f")

            open("$root/$f.dimacs.clean", "w") do io
                n = size(g.indexes, 1)

                for i = g.indexes
                    write(io, "$i\n")
                end

                for i = 1:size(g.indexes, 1)
                    u = g.indexes[i]
                    for j = i + 1:size(g.indexes, 1)
                        v = g.indexes[j]
                        if g.weights[u, v] > 0
                            write(io, "$u $v $(g.weights[u, v])\n")
                        end
                    end
                end
            end
        end
    end
end