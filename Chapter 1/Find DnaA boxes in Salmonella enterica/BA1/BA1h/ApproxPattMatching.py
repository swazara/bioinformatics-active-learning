#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1H
Title: Find All Approximate Occurrences of a Pattern in a String
URL: https://rosalind.info/problems/ba1h/

Description:
This script finds and returns all starting positions where the
input pattern appears with at most d mismatches.


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

def approximate_pattern_matching(pattern, dna_seq, d):
    d_match_ind = []
    patt_length = len(pattern)
    for i in range(len(dna_seq) - patt_length):
        h_distance = hamming_distance(pattern, dna_seq[i : i + patt_length])
        if h_distance <= d:
            d_match_ind.append(i)
    
    return d_match_ind



def main():
    if len(sys.argv) != 2:
        print("Use: python ApproxPattMatching.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            patt = lines[0].strip()
            seq = lines[1].strip()
            d = lines[2].strip()

            d_match_list = approximate_pattern_matching(patt, seq, int(d))
            print(" ".join(map(str, d_match_list)))
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()