#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3f
Title: Find an Eulerian Cycle in a Graph
URL: https://rosalind.info/problems/BA3f/

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

def find_eulerian_cycle(graph):
    temp_graph = {u: list(v) for u, v in graph.items()}

    start_node = list(temp_graph.keys())[0]
    stack = [start_node]
    cycle = []

    while stack:
        curr_node = stack[-1]
        if temp_graph.get(curr_node):
            next_node = temp_graph[curr_node].pop()
            stack.append(next_node)
        else: cycle.append(stack.pop())

    return cycle[::-1]

def main():
    if len(sys.argv) < 2:
        print("Usage: python FindEulerianCycle.py <input_file.txt>")
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

        eulerian_cycle = find_eulerian_cycle(parse_adjacency_list(lines))
        print("->".join([str(node) for node in eulerian_cycle]))


if __name__ == "__main__":
    main()