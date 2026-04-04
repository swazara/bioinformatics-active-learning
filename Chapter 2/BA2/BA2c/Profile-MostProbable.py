#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2c
Title: Find a Profile-most Probable k-mer in a String
URL: https://rosalind.info/problems/BA2c/

Description:
Given a DNA sequence and a profile matrix this algorithm finds the Profile-most
probable k-mer in the sequence.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
import numpy as np

def profile_most_probable_kmer(dna_seq, k, p_matrix):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    max_prob = -1.0
    most_prabable_kmer = ''
    
    for i in range(len(dna_seq) - k + 1):
        curr_prob = 1.0
        kmer = dna_seq[i:i+k]

        for j, base in enumerate(kmer):
            row = mapping[base]
            curr_prob *= p_matrix[row, j]
        
        if curr_prob > max_prob:
            max_prob = curr_prob
            most_prabable_kmer = kmer
    
    return most_prabable_kmer


def main():
    if len(sys.argv) != 2:
        print("Usage: python Profile-MostProbable.py <input_file.txt>")
        return

    file_path = sys.argv[1]
    data = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            
            data = (lines[0], # DNA Sequence
                    int(lines[1]), # k
                    np.array([line.split() for line in lines[2:]], dtype=float)) # Profile matrix
            
    except FileNotFoundError:
        print(f"Error: El archivo '{file_path}' no existe.")
    except Exception as e:
        print(f"Error al leer el archivo: {e}")

    if data:
        dna, k, p_matrix = data
        result = profile_most_probable_kmer(dna, k, p_matrix)
        print(result)

if __name__ == "__main__":
    main()