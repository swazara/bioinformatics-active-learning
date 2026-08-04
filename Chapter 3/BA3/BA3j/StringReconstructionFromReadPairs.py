#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3j
Title: Reconstruct a String from its Paired Composition
URL: https://rosalind.info/problems/BA3j/

Exception: This solution does not explore alternative Eulerian paths if the first one 
found is not valid (strings spelled dont overlap each other).

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
from collections import defaultdict

def RP_DeBrujin(kd_mers):
    # Constructs De Brujin graph from Read Pairs ((k,d)-mers) given as tuples (pattern_1, pattern_2)
    db_graph = defaultdict(list)
    for kdmer in kd_mers:
        db_graph[(kdmer[0][:-1],kdmer[1][:-1])].append((kdmer[0][1:],kdmer[1][1:]))
    
    return db_graph

def find_eulerian_path(graph):
    # Graph must be a digraph, strongly connected and have a solution 
    temp_graph = {u: list(v) for u, v in graph.items()}

    in_degrees = {}
    out_degrees = {}
    
    for source, targets in graph.items():
        out_degrees[source] = len(targets)
        for target in targets:
            in_degrees[target] = in_degrees.get(target, 0) + 1

    start_node = list(temp_graph.keys())[0]

    for node in out_degrees:
        # The starting node must have out_degrees(node) > in_degrees(node)
        if out_degrees[node] > in_degrees.get(node, 0):
            start_node = node
            break

    stack = [start_node]
    path = []

    while stack:
        curr_node = stack[-1]
        if temp_graph.get(curr_node):
            next_node = temp_graph[curr_node].pop()
            stack.append(next_node)
        else: path.append(stack.pop())
    
    return path[::-1]

def string_spelled_by_gapped_patterns(k,d,path):
    initial_node = path[0]
    string_prefix = initial_node[0]
    string_suffix = initial_node[1]

    for node in path[1:]:
        string_prefix += node[0][-1]
        string_suffix += node[1][-1]

    for i in range(k+d, len(string_prefix)):
        if string_prefix[i] != string_suffix[i-k-d]: return ""

    string_final = string_prefix + string_suffix[-(k + d):]
    return string_final

def main():
    if len(sys.argv) != 2:
        print("Usage: python StringReconstructionFromReadPairs.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    if lines:
        k, d = tuple(int(e) for e in lines[0].split())
        # Convert to tuples 
        kd_mers = [tuple(kdmer.split("|")) for kdmer in lines[1:]]

        db_graph = RP_DeBrujin(kd_mers)
        kdmers_path = find_eulerian_path(db_graph)
    
        print(string_spelled_by_gapped_patterns(k,d,kdmers_path))


if __name__ == "__main__":
    main()