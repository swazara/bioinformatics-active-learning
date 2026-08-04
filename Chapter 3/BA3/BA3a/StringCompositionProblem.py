#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3a
Title: Generate the k-mer Composition of a String
URL: https://rosalind.info/problems/BA3a/

Description:
Returns all the substrings of length k inside a string.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys


def generate_composition(string, k, order = True):
    l = len(string)
    composition = []

    for i in range(l - k + 1):
        composition.append(string[i:i+k])

    if order:
        return sorted(composition)
    else: return composition

    


def main():
    if len(sys.argv) < 2:
        print("Usage: python StringCompositionProblem.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
            data = (lines[0], # k
                    lines[1]) # DNA String

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    if data:
        k, dna_string = data
        comp = generate_composition(dna_string, int(k)) 
        
        for kmer in comp: print(kmer)

if __name__ == "__main__":
    main()