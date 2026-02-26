#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1E
Title: Find Patterns Forming Clumps in a Real Genome
URL: https://rosalind.info/problems/ba1e/

Description:
This script identifies (k, L, t)-clumps in a genomic sequence. A k-mer forms 
a clump if it appears at least 't' times within a window of length 'L'. 

Find regions with an unusually high frequency of specific patterns, 
such as DnaA boxes in the Origin of Replication.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
from collections import defaultdict

def get_frequency_table(text, k):
    #Creates a map of k-mers with the list of their starting positions.
    freq_map = defaultdict(list)
    for i in range(len(text) - k + 1):
        pattern = text[i : i + k]
        freq_map[pattern].append(i)
    return freq_map

def find_clumps(sequence, k, l_window, t):
    #Finds all k-mers that form a (k, L, t)-clump in the sequence.
    kmers_positions = get_frequency_table(sequence, k)
    clumps = []

    for kmer, positions in kmers_positions.items():
        if len(positions) >= t:
            for i in range(len(positions) - t + 1):
                distance = positions[i + t - 1] - positions[i]
                
                if distance <= l_window - k:
                    clumps.append(kmer)
                    break 

    return clumps

def main():
    if len(sys.argv) != 2:
        print("Usage: python ClumpFinder.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()    
            sequence = lines[0].strip()
            k, l_window, t = map(int, lines[1].split())
            
            result_clumps = find_clumps(sequence, k, l_window, t)

            print(" ".join(result_clumps))
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()