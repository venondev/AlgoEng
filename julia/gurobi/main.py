# Copyright 2020, Gurobi Optimization, LLC

# This example formulates and solves the following simple MIP model:
#  maximize
#        x +   y + 2 z
#  subject to
#        x + 2 y + 3 z <= 4
#        x +   y       >= 1
#        x, y, z binary

import sys
import math
import gurobipy as gp
from gurobipy import GRB

def main():

    # Create a new model
    m = gp.Model("mip1")

    n = int(sys.stdin.readline())
    cur = sys.stdin.readline()

    variables = {}
    has_edge = {}
    opt = None
    nodes = set()
    while cur != "\n":
        u, v, weight = cur.split(" ")
        u, v, weight = int(u), int(v), float(weight)

        isInfWeight = weight == -math.inf
        
        nodes.add(u)
        nodes.add(v)

        key = (u, v)
        var = m.addVar(
            lb=0, 
            ub=0 if isInfWeight else 1,
            vtype=GRB.CONTINUOUS,
            name=f"x_{0}_{1}"
        )

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

    # Edge Case if every edge has -inf edge weight
    if opt == None:
        print('Obj: 0')
        return

    nodes = list(nodes)
    nodes.sort()
    len_nodes = len(nodes)
    c_id = 0
    for u_idx in range(len_nodes):
        u = nodes[u_idx]
        for v_idx in range(u_idx + 1, len_nodes):
            v = nodes[v_idx]
            for w_idx in range(v_idx + 1, len_nodes):
                w = nodes[w_idx]

                c_id += 1
                m.addConstr(variables[(u, v)] +  variables[(v, w)] - variables[(u, w)] <= 1, str(c_id))

                c_id += 1
                m.addConstr(variables[(u, v)] -  variables[(v, w)] + variables[(u, w)] <= 1, str(c_id))

                c_id += 1
                m.addConstr(- variables[(u, v)] +  variables[(v, w)] + variables[(u, w)] <= 1, str(c_id))

    # for minimization
    # Set objective
    m.setObjective(opt, GRB.MINIMIZE)

    m.optimize()

    print('Obj: %g' % m.objVal)

main()