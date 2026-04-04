#!/usr/bin/env python3
"""
Rosalind Problem ID: BA2g
Title: Implement GibbsSampler
URL: https://rosalind.info/problems/BA2g/

Description:
Given a collection of DNA sequences, this algorithm starts with 
a random set of motifs (one per sequence). It then iteratively 
selects one sequence at random, removes its motif, and builds a 
profile matrix from the remaining motifs. From this profile, a 
new motif is sampled from the excluded sequence based on its 
probability distribution. This repeats N times, saving the best 
set of motifs that minimize the Score.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys
import numpy as np
import random

def profile_matrix(motifs):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    k = len(motifs[0])
    t = len(motifs)
    p_matrix = np.ones((4, k), dtype=float)

    for motif in motifs:
        for j in range(k):
            p_matrix[mapping[motif[j]]][j] += 1

    return p_matrix / (t + 4)

def score(motifs):
    k = len(motifs[0])
    t = len(motifs)
    total_score = 0
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    for j in range(k):
        counts = np.zeros(4)
        for i in range(t):
            counts[mapping[motifs[i][j]]] += 1
        total_score += (t - np.max(counts))
    return int(total_score)

def get_prob(kmer, p_matrix):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    prob = 1.0
    for j, base in enumerate(kmer):
        prob *= p_matrix[mapping[base]][j]
    return prob

def profile_random_kmer(dna_i, k, p_matrix):
    n = len(dna_i)
    kmers = [dna_i[j:j+k] for j in range(n - k + 1)]
    probs = np.array([get_prob(km, p_matrix) for km in kmers])
    probs /= probs.sum()
    idx = np.random.choice(len(kmers), p=probs)
    return kmers[idx]

def gibbs_sampler(dna, k, t, N):
    motifs = []
    for seq in dna:
        r = random.randint(0, len(dna[0]) - k)
        motifs.append(seq[r:r+k])
    
    best_motifs = list(motifs)
    
    for _ in range(N):
        i = random.randint(0, t - 1)
        reduced_motifs = [motifs[j] for j in range(t) if j != i]
        p_matrix = profile_matrix(reduced_motifs)
        motifs[i] = profile_random_kmer(dna[i], k, p_matrix)
        
        if score(motifs) < score(best_motifs):
            best_motifs = list(motifs)
            
    return best_motifs, score(best_motifs)

def main():
    if len(sys.argv) != 2:
        print("Usage: python GibbsSampler.py <input_file.txt>")
        return
    file_path = sys.argv[1]

    try:
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]
            k, t, N = map(int, lines[0].split())
            dna_strings = lines[1:]
        
        best_overall = None
        min_score = float('inf')
        
        for _ in range(20):
            current_motifs, current_score = gibbs_sampler(dna_strings, k, t, N)
            if current_score < min_score:
                min_score = current_score
                best_overall = current_motifs
        
        print("\n".join(best_overall))
        
    except FileNotFoundError:
        print(f"Error: El archivo '{file_path}' no existe.")
    except Exception as e:
        print(f"Error al leer el archivo: {e}")


if __name__ == "__main__":
    main()