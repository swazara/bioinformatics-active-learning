#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1n
Title: Generate the d-Neighborhood of a String
URL: https://rosalind.info/problems/ba1n/

Description:
This script returns all the k-mers within at most a *d* hamming distance 
from a given pattern (of length k). This is also known as the d-Neighborhood of a pattern.


Author: Santiago Wilders Azara 
Date: 2026
"""

import sys
def d_neighborhood(pattern, d):

    def generate(idx, remaining_d):
        if idx == len(pattern):
            return [""]
        results = []
        original = pattern[idx]
        for suffix in generate(idx + 1, remaining_d):
            results.append(original + suffix)
        if remaining_d > 0:
            for base in "ACGT":
                if base != original:
                    for suffix in generate(idx + 1, remaining_d - 1):
                        results.append(base + suffix)
        return results

    return set(generate(0, d))



def main():
    if len(sys.argv) != 2:
        print("Use: python GenerateNeighborhood.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            patt = lines[0].strip()
            d = lines[1].strip()

            d_neigh = d_neighborhood(patt, int(d))
            print("\n".join(map(str, d_neigh)))
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()