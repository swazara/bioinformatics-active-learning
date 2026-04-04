#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2h
Title: Implement DistanceBetweenPatternAndStrings
URL: https://rosalind.info/problems/BA2h/

Description:
Returns the minimum accumulated distances from a pattern to a set of DNA sequences.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys

def hamming_distance(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))


def distance_between_pattern_and_strings(pattern, dna_sequences):
    k = len(pattern)
    distance = 0
    for seq in dna_sequences:
            min_hd = float('inf')
            for j in range(len(seq) - k + 1):
                hd = hamming_distance(pattern, seq[j: j+k])
                if hd < min_hd: min_hd = hd
            distance += min_hd
    return distance

def main():
    if len(sys.argv) < 2:
        print("Usage: python DistanceBetweenPatternAndStrings.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
            data = (lines[0], # pattern
                    lines[1].split()) # DNA String

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    if data:
        pattern, dna_sequences = data
        dist = distance_between_pattern_and_strings(pattern, dna_sequences) 
        print(dist)

if __name__ == "__main__":
    main()