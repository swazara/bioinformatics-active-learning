#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3e
Title: Construct the De Bruijn Graph of a Collection of k-mers
URL: https://rosalind.info/problems/BA3e/

Description:
Returns a De Brujin graph constructed from a collection of k-mers


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

def main():
    if len(sys.argv) < 2:
        print("Usage: python DeBrujinPatterns.py <input_file.txt>")
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
        db_graph = DeBrujin(lines)

        for source, targets in sorted(db_graph.items()):
            targets_string = ",".join(sorted(targets))
            print(source + " -> " + targets_string)
            


if __name__ == "__main__":
    main()