#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3c
Title: Construct the Overlap Graph of a Collection of k-mers
URL: https://rosalind.info/problems/BA3c/

Description:
Constructs an overlap graph for a given set of patterns. Each pattern/k-mer is a node,
if first k − 1 nucleotides of a k-mer and last k − 1 nucleotides of another k-mer are equal, 
an edge is defined between them.


Author: Santiago Wilders Azara
Date: 2026
"""

import sys


def prefix(pattern):
    return pattern[:-1]
    
def suffix(pattern):
    return pattern[1:]

def overlap_graph(patterns):
    edge_list = [] # It's an edge list even though Rosalind considers it an adjacency list.

    for i in range(len(patterns)):
        p_i = patterns[i]
        for j in range(i+1, len(patterns)):
            p_j = patterns[j]

            if prefix(p_i) == suffix(p_j):   
                edge_list.append((p_j, p_i))

            if suffix(p_i) == prefix(p_j):
                edge_list.append((p_i, p_j))

    return edge_list


def main():
    if len(sys.argv) < 2:
        print("Usage: python OverlapGraph.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            patterns = [line.strip() for line in f.readlines() if line.strip()]

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    if patterns:
        o_graph = overlap_graph(patterns)
        for edge in o_graph:
            print(edge[0] + " -> " + edge[1])

if __name__ == "__main__":
    main()