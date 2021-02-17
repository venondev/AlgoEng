using JuMP, Suppressor, Gurobi

include("../helper/graph.jl")
include("../helper/p3.jl")
include("../data_reduction/data_reduction.jl")

const GRB_ENV = @suppress Gurobi.Env()

function solveLP(g::Graph)
    model = JuMP.direct_model(Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "OutputFlag", 1)

    numVars = 0
    for i = 1:size(g.indexes, 1)
        for j = i + 1:size(g.indexes, 1)
            numVars += 1
        end
    end
    @variable(model, 1 >= x[1:numVars] >= 0)

    vars = Dict{Tuple{Int, Int}, Int}()
    obj = nothing
    edgeIdx = 1
    for i = 1:size(g.indexes, 1)
        u = g.indexes[i]
        for j = i + 1:size(g.indexes, 1)
            v = g.indexes[j]

            vars[(u, v)] = edgeIdx

            w = abs(getWeight(g, u, v))

            if w == Inf
                fix(x[edgeIdx], 0., force=true)
                edgeIdx += 1
                continue
            end

            if hasEdge(g, u, v)
                t = w * (1 - x[edgeIdx])
                if obj === nothing
                    obj = t
                else
                    obj += t
                end
            else 
                t2 = w * x[edgeIdx]
                if obj === nothing
                    obj = t2
                else
                    obj += t2
                end
            end

            edgeIdx += 1
        end
    end

    if obj === nothing
        return -1 # TODO: Nochmal überdenken vllt sogar Inf? 
    end

    for i = 1:size(g.indexes, 1)
        u = g.indexes[i]
        for j = i + 1:size(g.indexes, 1)
            v = g.indexes[j]
            for k = j + 1:size(g.indexes, 1)
                w = g.indexes[k]

                a = x[vars[(u, v)]]
                b = x[vars[(v, w)]]
                c = x[vars[(u, w)]]

                @constraint(model, a + b - c <= 1.)
                @constraint(model, a - b + c <= 1.)
                @constraint(model,-a + b + c <= 1.)
            end
        end
    end
    
    @objective(model, Min, obj)

    optimize!(model)
    
    return objective_value(model)
end