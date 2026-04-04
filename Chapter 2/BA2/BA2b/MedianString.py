#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2b
Title: Implement MotifEnumeration
URL: https://rosalind.info/problems/BA2b/

Description:
Finds a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers 'Pattern', 
where d(Pattern, Dna) is the minimum hamming distance between 'Pattern' and an k-mer in 'Dna',
a set of 't' sequences with 'n' nucleotides of length.  We call such a k-mer a median string for Dna.
This is a pretty expensive algorithm with O(4^k.n.k.t) complexity.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys

def hamming_distance(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def numberToPattern(index, k):
    bases = ['A', 'C', 'G', 'T']
    pattern = ""
    for _ in range(k):
        pattern = bases[index % 4] + pattern
        index //= 4
    return pattern

def median_string(k, dna_sequences):
    min_count = float('inf')
    median = ''
    
    for i in range(4**k):
        pattern_i = numberToPattern(i, k)
        curr_count = 0
        
        for seq in dna_sequences:
            min_hd = float('inf')
            for j in range(len(seq) - k + 1):
                hd = hamming_distance(pattern_i, seq[j: j+k])
                if hd < min_hd: min_hd = hd
            curr_count += min_hd
        
        if curr_count < min_count:
            min_count = curr_count
            median = pattern_i
    
    return median


def main():
    if len(sys.argv) < 2:
        print("Usage: python MedianString.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
            data = (lines[0], # k
                    lines[1:len(lines)]) # DNA Strings

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    if data:
        k, dna_strings = data
        median = median_string(int(k), dna_strings)
        print(median)

if __name__ == "__main__":
    main()