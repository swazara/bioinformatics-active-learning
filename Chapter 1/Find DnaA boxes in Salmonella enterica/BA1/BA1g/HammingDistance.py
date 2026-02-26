#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1G
Title: Compute the Hamming Distance Between Two Strings
URL: https://rosalind.info/problems/ba1g/

Description:
This script finds mismatches between two strings and then returns 
the total, also known as the Hamming Distance.

In bioinformatics, this metric is essential for identifying SNPs 
(Single Nucleotide Polymorphisms) and quantifying point mutations 
between homologous sequences.

Author: Santiago Wilders Azara 
Date: 2026
"""

import sys

def hamming_distance(p,q):
    count = 0
    for i in range(len(p)):
        if p[i] != q [i]:
            count+=1

    return count

def main():
    if len(sys.argv) != 2:
        print("Use: python HammingDistance.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            p = lines[0].strip()
            q = lines[1].strip()

            if len(p) != len(q):
                return
            
            print(hamming_distance(p,q))
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()