#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1F
Title: Find All Local Minima of the Skew Diagram
URL: https://rosalind.info/problems/ba1f/

Description:
This script identifies all positions 'i' (from 0 to |Genome|) where the 
Skew_i(Genome) achieves its minimum value. The skew is the cumulative 
difference between Guanine (G) and Cytosine (C), #G - #C.

In bacterial genomics, the global minimum of the skew diagram is an 
approach to predict the Origin of Replication (OriC).

Author: Santiago Wilders Azara 
Date: 2026
"""

import sys

def minimumSkew(dna_seq):
    skew_values = {'G': 1, 'C': -1, 'A': 0, 'T': 0}
    curr_val = 0
    min_val = 0
    ind_list = [0] # The sequence may be starting at a minimum
    
    for i, base in enumerate(dna_seq.upper(), 1):
        curr_val += skew_values.get(base, 0)
        
        if curr_val < min_val:
            min_val = curr_val
            ind_list = [i]
        elif curr_val == min_val:
            ind_list.append(i)
            
    return ind_list

def main():
    if len(sys.argv) != 2:
        print("Use: python script.py file.txt")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            dna_seq = f.read().strip()
            
            skew_mins = minimumSkew(dna_seq)

            print(" ".join(map(str, skew_mins)))
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()