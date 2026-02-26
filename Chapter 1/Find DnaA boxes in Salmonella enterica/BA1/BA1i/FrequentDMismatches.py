#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1i
Title: Find the Most Frequent Words with Mismatches in a String
URL: https://rosalind.info/problems/ba1i/

Description:
This script receives a string (dna sequence) as well as integers k and d,
then finds and returns the most frequent k-mers with up to d mismatches 
in the sequence.

Author: Santiago Wilders Azara 
Date: 2026
"""

import sys

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
    

def frequent_pattern_mismatches(dna_seq, k, d):
    frequency_map = {}
    most_freq = 0
    most_freq_patt = set()
    for i in range(len(dna_seq) - k + 1):
        pattern = dna_seq[i : i+k]
        neighbors = neighborhood(pattern,d)
        for n in neighbors:
            frequency_map[n] = frequency_map.get(n, 0) + 1
            if  frequency_map[n] > most_freq:
                most_freq_patt = {n}
                most_freq = frequency_map[n]

            elif frequency_map[n] == most_freq:
                most_freq_patt.add(n)
    
    #print(sorted(frequency_map.items(), key= lambda item: item[1],reverse=True))
    return most_freq_patt

def main():
    if len(sys.argv) != 2:
        print("Use: python FrequentDMismatches.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            dna_seq = lines[0].strip()
            k, d = (lines[1].strip()).split()
            d_match_set = frequent_pattern_mismatches(dna_seq, int(k), int(d))
            
            print(" ".join(map(str, d_match_set)))
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()