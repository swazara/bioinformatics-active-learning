#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3h
Title: Reconstruct a String from its k-mer Composition
URL: https://rosalind.info/problems/BA3h/

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
from collections import defaultdict

def DeBrujin(patterns):
    db_graph = defaultdict(list)
    for kmer in patterns:
        db_graph[kmer[:-1]].append(kmer[1:])
    
    return db_graph

def find_eulerian_path(graph):
    # graph must be a digraph, strongly connected and have a solution 
    temp_graph = {u: list(v) for u, v in graph.items()}

    in_degrees = {}
    out_degrees = {}
    
    for source, targets in graph.items():
        out_degrees[source] = len(targets)
        for target in targets:
            in_degrees[target] = in_degrees.get(target, 0) + 1

    start_node = list(temp_graph.keys())[0]

    for node in out_degrees:
        # The starting node must be the one with more out than in edges.
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

def string_spelled_from_kmers(kmers):
    if not kmers:
        return
    
    n = len(kmers) 
    k = len(kmers[0])
    string_spelled = kmers[0]

    for i in range(1, n):
        string_spelled += kmers[i][-1]

    return string_spelled

def main():
    if len(sys.argv) < 2:
        print("Usage: python StringReconstruction.py <input_file.txt>")
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
        k = lines[0]
        patterns = lines[1:]

        db_graph = DeBrujin(patterns)
        kmers_path = find_eulerian_path(db_graph)

        print(string_spelled_from_kmers(kmers_path))


if __name__ == "__main__":
    main()