#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3m
Title: Generate All Maximal Non-Branching Paths in a Graph
URL: https://rosalind.info/problems/BA3m/

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
from collections import defaultdict

def parse_adjacency_list(lines):
    graph = defaultdict(list)
        
    for line in lines:
        source_str, targets_str = line.split(" -> ")

        source = int(source_str)
        targets = [int(t) for t in targets_str.split(",")]
        
        graph[source].extend(targets)
            
    return dict(graph)

def maximal_non_branching_paths(graph):
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    all_nodes = set()

    for node, edges in graph.items():
        all_nodes.add(node)
        out_degree[node] = len(edges)
        for node_out in edges:
            in_degree[node_out] += 1
            all_nodes.add(node_out)

    paths = []
    visited = set()

    for node in all_nodes:
        if in_degree[node] != 1 or out_degree[node] != 1:
            if out_degree[node] > 0:
                for node_out in graph[node]:
                    non_branching_path = [node, node_out]
                    current = node_out
                    
                    while in_degree[current] == 1 and out_degree[current] == 1:
                        visited.add(current)
                        next_node = graph[current][0]
                        non_branching_path.append(next_node)
                        current = next_node
                        
                    paths.append(non_branching_path)

    for node in all_nodes:
        if in_degree[node] == 1 and out_degree[node] == 1 and node not in visited:
            visited.add(node)
            cycle = [node]
            current = graph[node][0]
            
            while current != node:
                visited.add(current)
                cycle.append(current)
                current = graph[current][0]
                
            cycle.append(node)
            paths.append(cycle)

    return paths
                        


def main():
    if len(sys.argv) != 2:
        print("Usage: python MaximalNonBranchingPaths.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return
    except Exception as e:
        print(f"Unexpected error: {e}")
        return
    
    input_graph = parse_adjacency_list(lines)
    non_batching_paths = maximal_non_branching_paths(input_graph)

    for path in non_batching_paths:
        print(" -> ".join([str(node) for node in path]))


if __name__ == "__main__":
    main()