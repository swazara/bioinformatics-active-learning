#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2f
Title: Implement RandomizedMotifSearch
URL: https://rosalind.info/problems/BA2f/

Description:
This is a Monte Carlo algorithm, a randomized approach that does not guarantee an 
optimal solution for the Subtle Motif Problem in a single run. Given a collection 
of DNA sequences, the algorithm starts from a random set of motifs (one from each 
sequence). It then iteratively improves this set to find the best motifs that 
minimize the Score (the sum of differences within motifs).

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
import numpy as np
import random

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
    p_matrix = np.ones((4,k),dtype=float)

    for i in range(n):
        for j in range(k):
            p_matrix[mapping[motifs[i][j]]][j] += 1

    return p_matrix / (n + 4)
     
        
def score(motifs_selection):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    n = len(motifs_selection)
    k = len(motifs_selection[0])
    total_score = 0

    for j in range(k):
        base_count = np.zeros(4, dtype=int)
        for i in range(n):
            base_count[mapping[motifs_selection[i][j]]] += 1
        
        max_idx = 0
        for b in range(1,4):
            if base_count[b] > base_count[max_idx]:
                max_idx = b

        total_score += n - base_count[max_idx]

    return total_score

def motifs(p_matrix, dna):
    k = len(p_matrix[0])
    motif_selection = np.array(dna, dtype=f"U{k}")
    t = len(dna)

    for i in range(t):
        seq_kmer = profile_most_probable_kmer(dna[i], k, p_matrix)
        motif_selection[i] = seq_kmer
    
    return motif_selection

def randomized_motif_search(dna, k, t):
    motif_selection = np.array(dna, dtype=f"U{k}")
    t = len(dna)
    n = len(dna[0])

    for i in range(t):
        random_idx = random.randint(0, n-k)
        motif_selection[i] = dna[i][random_idx:random_idx+k]

    best_motifs = motif_selection.copy()
    
    while True:
        profile = profile_matrix(motif_selection)
        motif_selection = motifs(profile, dna)
        
        if score(motif_selection) < score(best_motifs):
            best_motifs = motif_selection.copy()
        else:
            return best_motifs
 
def main():
    if len(sys.argv) != 2:
        print("Usage: python RandomizedMotifSearch.py <input_file.txt>")
        return

    file_path = sys.argv[1]
    data = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            k, t =  lines[0].split()

            data = (int(k),
                    int(t),
                    lines[1:])
            
    except FileNotFoundError:
        print(f"Error: El archivo '{file_path}' no existe.")
    except Exception as e:
        print(f"Error al leer el archivo: {e}")

    if data:
        k, t, dna_strings = data
        
        best_motifs = randomized_motif_search(dna_strings, k, t)
        best_score = score(best_motifs)
        
        for _ in range(999):
            current_motifs = randomized_motif_search(dna_strings, k, t)
            current_score = score(current_motifs)
            
            if current_score < best_score:
                best_score = current_score
                best_motifs = current_motifs

        print("\n".join(best_motifs))

if __name__ == "__main__":
    main()