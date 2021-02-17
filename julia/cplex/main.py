import docplex.mp.model as cpx
import sys
import math


def main():
    opt_model = cpx.Model(name="MIP Model")
    opt_model.context.cplex_parameters.workmem=1024

    n = int(sys.stdin.readline())

    cur = sys.stdin.readline()

    variables = {}
    has_edge = {}
    opt = None
    nodes = set()
    while cur != "":
        u, v, weight = cur.split(" ")
        u, v, weight = int(u), int(v), float(weight)

        isInfWeight = weight == -math.inf
        
        nodes.add(u)
        nodes.add(v)

        key = (u, v)
        var = opt_model.integer_var(lb=0, ub=1, name="x_{0}_{1}".format(u,v))

        if not isInfWeight:
            weight = int(weight)
            cur_opt = abs(weight) * var
            if weight > 0:
                cur_opt = (1 - var) * weight

            if opt is None:
                opt = cur_opt
            else:
                opt += cur_opt

        variables[key] = var
        has_edge[key] = weight > 0
        cur = sys.stdin.readline()

    nodes = list(nodes)
    len_nodes = len(nodes)
    c_id = 0
    for u_idx in range(len_nodes):
        u = nodes[u_idx]
        for v_idx in range(u_idx + 1, len_nodes):
            v = nodes[v_idx]
            for w_idx in range(v_idx + 1, len_nodes):
                w = nodes[w_idx]

                c_id += 1
                opt_model.add_constraint(ct=(variables[(u, v)] +  variables[(v, w)] - variables[(u, w)] <= 1), ctname=str(c_id))

                c_id += 1
                opt_model.add_constraint(ct=(variables[(u, v)] -  variables[(v, w)] + variables[(u, w)] <= 1), ctname=str(c_id))

                c_id += 1
                opt_model.add_constraint(ct=(- variables[(u, v)] +  variables[(v, w)] + variables[(u, w)] <= 1), ctname=str(c_id))

    # for minimization
    opt_model.minimize(opt)
    has_sol = opt_model.solve()
    

    if has_sol:
        for key in variables:
            sol_has_edge = int(float(variables[key])) == 1
            if (sol_has_edge and not has_edge[key]) or (not sol_has_edge and has_edge[key]):
                print(str(key[0]) + " " + str(key[1]))

main()