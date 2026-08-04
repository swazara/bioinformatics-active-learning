#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3d
Title: Construct the De Bruijn Graph of a String
URL: https://rosalind.info/problems/BA3d/

Description:
Returns a De Brujin graph represented by a dictionary of lists. The edges soted before printing.


Author: Santiago Wilders Azara
Date: 2026
"""

import sys
from collections import defaultdict

def DeBrujin(dna_string, c = 3):
    db_graph = defaultdict(list)
    k = c - 1

    for i in range(len(dna_string) - k):
        j = i +1
        db_graph[dna_string[i:i+k]].append(dna_string[j:j+k])
    
    return db_graph

def main():
    if len(sys.argv) < 2:
        print("Usage: python DeBrujinGraph.py <input_file.txt>")
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
        k = int(lines[0])
        dna_string = lines[1]

        db_graph = DeBrujin(dna_string, k)

        for source, targets in sorted(db_graph.items()):
            targets_string = ", ".join(sorted(targets))
            print(source + " -> " + targets_string)
            


if __name__ == "__main__":
    main()