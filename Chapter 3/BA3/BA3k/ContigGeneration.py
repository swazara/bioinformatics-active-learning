#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3k
Title: Generate Contigs from a Collection of Reads
URL: https://rosalind.info/problems/BA3k/

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

def reconstruct_contig(contig_path):
    # contig_path must be a string list
    contig = contig_path[0]
    for i in range(1,len(contig_path)):
        contig += contig_path[i][-1]

    return contig



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
        print("Usage: python ContigGeneration.py <input_file.txt>")
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
    
    db_graph = DeBrujin(lines)
    non_batching_paths = maximal_non_branching_paths(db_graph)

    res = map(reconstruct_contig, non_batching_paths)
    
    print(" ".join(res))


if __name__ == "__main__":
    main()