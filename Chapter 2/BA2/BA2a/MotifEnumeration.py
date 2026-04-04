#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2a
Title: Implement MotifEnumeration
URL: https://rosalind.info/problems/BA2a/

Description:
Searches for (k,d)-motifs in a collection of DNA strings. With 'k' being the k-mer length with at most 'd' mismatches.
This is a computationally expensive brute force algorithm.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
import numpy as np

def hamming_distance(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def get_val(symbol):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return mapping.get(symbol.upper(), 0)

def patternToNumber(pattern):
    number = 0
    for char in pattern:
        number = (number * 4) + get_val(char)
    return number

def numberToPattern(index, k):
    bases = ['A', 'C', 'G', 'T']
    pattern = ""
    for _ in range(k):
        pattern = bases[index % 4] + pattern
        index //= 4
    return pattern

def neighborhood(pattern, d):

    def generate(current_idx, current_d):
        if current_idx == len(pattern):
            return [""]
        
        results = []
        char_original = pattern[current_idx]
        
        for suffix in generate(current_idx + 1, current_d):
            results.append(char_original + suffix)
            
        if current_d > 0:
            for char in ["A", "C", "G", "T"]:
                if char != char_original:
                    for suffix in generate(current_idx + 1, current_d - 1):
                        results.append(char + suffix)
        return results

    return set(generate(0, d))

def motif_enumeration(k, d, dna_sequences):
    freqs = np.zeros(4**k)
    patterns = set()
    for i in range(4**k):
        pattern_i = numberToPattern(i, k)
        neigh_set = neighborhood(pattern_i, d)

        for seq in dna_sequences:
            for pattern_j in neigh_set:
                if pattern_j in seq:
                    freqs[i] += 1
                    break
        
        if freqs[i] == len(dna_sequences):
            patterns.add(pattern_i)

    return patterns



def main():
    if len(sys.argv) < 2:
        print("Usage: python MotifEnumeration.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            # BA2a input format: Line 0 = k d, Line 1 to End = DNA string
            k, d = lines[0].strip().split()
            dna_strings = lines[1:len(lines)]

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")

    patterns = motif_enumeration(int(k), int(d), dna_strings)
    print(" ".join(map(str, patterns)))
    
if __name__ == "__main__":
    main()