#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2d
Title: Implement GreedyMotifSearch
URL: https://rosalind.info/problems/BA2d/

Description:
Given a collection of DNA sequences, this algorithm returns the best 
motifs found based on basic score (sum of differences with consensus) 
and greedy decisions. This algorithm does not guarantee optimal solution 
for the Subtle Motif Problem.


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

def profile_matrix(motifs):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    k = len(motifs[0])
    n = len(motifs)
    p_matrix = np.zeros((4,k),dtype=float)

    for i in range(n):
        for j in range(k):
            p_matrix[mapping[motifs[i][j]]][j] += 1

    return p_matrix/n
     
        
def score(motifs):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    n = len(motifs)
    k = len(motifs[0])
    total_score = 0

    for j in range(k):
        base_count = np.zeros(4, dtype=int)
        for i in range(n):
            base_count[mapping[motifs[i][j]]] += 1
        
        max_idx = 0
        for b in range(1,4):
            if base_count[b] > base_count[max_idx]:
                max_idx = b

        total_score += n - base_count[max_idx]

    return total_score


def greedy_motif_search(k, t, dna):
    # Automatically getting the first k bases, Uk means Unicode data format of length k.
    best_motifs = np.array(dna, dtype=f"U{k}")

    for w in range(len(dna[0]) - k + 1):
        motifs = [dna[0][w:w+k]]
        for i in range(1,t):
            profile = profile_matrix(motifs)
            motifs.append(profile_most_probable_kmer(dna[i],k,profile))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs
            
    return best_motifs

def main():
    if len(sys.argv) != 2:
        print("Usage: python GreedyMotifSearch.py <input_file.txt>")
        return

    file_path = sys.argv[1]
    data = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            k, t =  lines[0].split()

            data = (int(k),
                    int(t),
                    lines[1:]) # Collection of strings
            
    except FileNotFoundError:
        print(f"Error: El archivo '{file_path}' no existe.")
    except Exception as e:
        print(f"Error al leer el archivo: {e}")

    if data:
        k, t, dna_strings = data
        result = greedy_motif_search(k, t, dna_strings)

        print("\n".join(result))

if __name__ == "__main__":
    main()