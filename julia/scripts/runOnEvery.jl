include("../main.jl")

path = "/home/hjal/CLionProjects/AlgoEng/students"

for (root, dirs, files) in walkdir(path)
    println("Checking $root")
    for f = files
        if endswith(f, ".dimacs")
            matrix = readInputFromFile("$root/$f")
            compareP3(matrix, f)
        end
    end
end