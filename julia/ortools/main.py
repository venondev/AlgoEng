from __future__ import print_function
from ortools.sat.python import cp_model

import sys
import fileinput
import math

def main():
    model = cp_model.CpModel()

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

        key = f"{u} {v}"
        var = model.NewIntVar(0, 0 if isInfWeight else 1, key)

        if not isInfWeight:
            weight = int(weight)
            cur_opt = abs(weight) * var
            if weight > 0:
                cur_opt = (1 - var) * weight

            if opt == None:
                opt = cur_opt
            else:
                opt += cur_opt

        variables[key] = var
        has_edge[key] = weight > 0
        cur = sys.stdin.readline()

    nodes = list(nodes)
    len_nodes = len(nodes)
    for u_idx in range(len_nodes):
        u = nodes[u_idx]
        for v_idx in range(u_idx + 1, len_nodes):
            v = nodes[v_idx]
            for w_idx in range(v_idx + 1, len_nodes):
                w = nodes[w_idx]

                model.Add(variables[f"{u} {v}"] + variables[f"{v} {w}"] - variables[f"{u} {w}"] <= 1)
                model.Add(variables[f"{u} {v}"] - variables[f"{v} {w}"] + variables[f"{u} {w}"] <= 1)
                model.Add(- variables[f"{u} {v}"] + variables[f"{v} {w}"] + variables[f"{u} {w}"] <= 1)

    
    model.Minimize(opt)

    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    if status == cp_model.OPTIMAL:
        for key in variables:
            print(key + " " + str(solver.Value(variables[key])))

if __name__ == '__main__':
  main()