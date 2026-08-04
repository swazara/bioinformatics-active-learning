#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3i
Title: k-Universal Circular String Problem
URL: https://rosalind.info/problems/BA3i/

Description: Finds the k-universal circular string constructing a De Brujin graph with 
(k-1)-mers as nodes and finds an Eulerian path.

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

def get_cycle_from_kmers(kmers):
    if not kmers:
        return
    
    n = len(kmers) 
    k = len(kmers[0])
    string_spelled = kmers[0]

    for i in range(1, n):
        string_spelled += kmers[i][-1]

    return string_spelled[:-k]

def main():
    if len(sys.argv) != 2:
        print("Usage: python kUniversalString.py [k]")
        return

    k = int(sys.argv[1])
    
    binary_kmers = [f"{i:0{k}b}" for i in range(2**k)]
    db_graph = DeBrujin(binary_kmers)
    kmers_path = find_eulerian_cycle(db_graph)

    print(get_cycle_from_kmers(kmers_path))


if __name__ == "__main__":
    main()