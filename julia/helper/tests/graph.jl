using Test
include("../graph.jl")

function smallMerge()
    g = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs")
    g_merged = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs.merged")

    merge!(g, 2, 3)
    g_merged.indexes = [1, 2, 4]

    return graphEqual(g, g_merged)
end

function smallUnmerge()
    g = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs")
    g_copy = copy(g)

    merge!(g, 2, 3)
    unmerge!(g, false)

    return graphEqual(g, g_copy)
end

function smallUnmergeChecked()
    g = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs")
    g_copy = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs")

    weight = getWeight(g, 1, 2)
    setWeight(g, 1, 2, weight > 0 ? Inf : -Inf)

    weight = getWeight(g, 1, 2)
    setWeight(g_copy, 1, 2, weight > 0 ? Inf : -Inf)

    k_decr = merge!(g, 2, 3)
    unmerge!(g, false)

    if k_decr != 6 return 
        println("Wrong k_decr")
        return false
    end

    return graphEqual(g, g_copy)
end

# Wie smallUnmergeChecked nur das 2 in 3 gemerged wird
function smallUnmergeCheckedMergeIn3()
    g = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs")
    g_copy = readMatrixFromFile("$(@__DIR__)/test_graphs/graph1.dimacs")

    weight = getWeight(g, 1, 2)
    setWeight(g, 1, 2, weight > 0 ? Inf : -Inf)

    weight = getWeight(g, 1, 2)
    setWeight(g_copy, 1, 2, weight > 0 ? Inf : -Inf)

    k_decr = merge!(g, 3, 2)
    unmerge!(g, false)

    if k_decr != 6 return 
        println("Wrong k_decr")
        return false
    end

    return graphEqual(g, g_copy)
end

function smallClusterRemove()
    g = readMatrixFromFile("$(@__DIR__)/test_graphs/cluster_graph1.dimacs")
    g_copy = readMatrixFromFile("$(@__DIR__)/test_graphs/cluster_graph1.dimacs")

    removed_nodes = removeCluster(g, [2])

    if sort(g.indexes) != [4, 5]
        return false
    end

    undoRemoveCluster(g, removed_nodes)

    return graphEqual(g, g_copy)
end

@testset "graph merge" begin
    @test smallMerge()
    @test smallUnmerge()
    @test smallUnmergeChecked()
    @test smallUnmergeCheckedMergeIn3()
    @test smallClusterRemove()
end