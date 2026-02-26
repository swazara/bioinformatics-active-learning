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

def hamming_distance(p,q):
    c = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            c += 1
    return c

def d_neighborhood(pattern, d):
    l = len(pattern)
    if d == 0:
        return set()
    if l == 1:
        return set({'A', 'C', 'G', 'T'})
    
    neighborhood = set()
    suffix_neighborhood = d_neighborhood(pattern[1:l], d)
    for s in suffix_neighborhood:
        if hamming_distance(s, pattern[1:l]) < d:
            for n in {'A', 'C', 'G', 'T'}:
                neighborhood.add(n+s)
        else:
            neighborhood.add(pattern[0]+s)
    return neighborhood



def main():
    if len(sys.argv) != 2:
        print("Use: python bookGenerateNeighborhood.py <input.txt>")
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