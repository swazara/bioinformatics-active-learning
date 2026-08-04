#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3l
Title: Construct a String Spelled by a Gapped Genome Path
URL: https://rosalind.info/problems/BA3l/

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
from collections import defaultdict

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
        print("Usage: python GappedGenomePathString.py <input_file.txt>")
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
        print(string_spelled_by_gapped_patterns(k,d,kd_mers))


if __name__ == "__main__":
    main()