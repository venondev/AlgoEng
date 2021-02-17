include("../helper/graph.jl")
include("merge_greater_k.jl")
include("heavy_edge.jl")
include("neighborhood.jl")
include("cluster_remove.jl")

function doDataReductionPre(g::Graph)::Float64
    removedClusters = removeClusters(g)
    N, V = findAdjWeights(g)
    k_decr = 0.
    merges = 0

    heavyNonEdge!(g, N)

    merges_temp::Int, k_decr_temp::Float64 = heavyEdgeSingleEnd!(g, 0., N, V)
    k_decr += k_decr_temp
    merges += merges_temp

    merges_temp_2::Int, k_decr_temp_2::Float64 = heavyEdgeBothEnds!(g, 0., N, V)
    k_decr += k_decr_temp_2
    merges += merges_temp_2

    merges_temp_3::Int, k_decr_temp_3::Float64 = large_neighborhood(g, 0., N, V)
    k_decr += k_decr_temp_3
    merges += merges_temp_3

    if merges > 0
        k_decr_temp = doDataReductionPre(g)
        k_decr += k_decr_temp
    end

    return k_decr
end

function doDataReductionBranch(g::Graph, k::Float64, recDepth::Int)::Tuple{Int,Int}
    if recDepth % 3 == 0 
        removedClusters = 0 #removeClusters(g)

        k_decr = 0.
        merges = 0

        merges_temp, k_decr_temp = mergeGreaterK(g, k)
        k_decr += k_decr_temp
        merges += merges_temp

        if merges > 0
            merges_temp, k_decr_temp = doDataReductionBranch(g, k - k_decr, recDepth)
            k_decr += k_decr_temp
            merges += merges_temp
        end

        return merges + removedClusters, k_decr
    end

    return (0, 0)
end

