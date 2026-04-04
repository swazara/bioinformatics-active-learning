#!/usr/bin/env python3
"""
motif_utils.py — Motif finding algorithms for DNA sequence analysis.

Unified motif-finding module consolidating the Chapter 2 algorithms from 
Bioinformatics Algorithms: An Active Learning Approach, with shared infrastructure 
and interchangeable scoring functions. Consolidation made with assistance from Claude.

Author: Santiago Wilders Azara
Date: 2026

Notation used in complexity comments
--------------------------------------
  n  : length of each DNA sequence
  k  : motif length
  t  : number of DNA sequences
  N  : number of Gibbs sampler iterations
"""

import sys
import numpy as np
import random

# ---------------------------------------------------------------------------
# Shared mapping
# ---------------------------------------------------------------------------

MAPPING = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# ---------------------------------------------------------------------------
# Scoring functions
# Two options available for greedy_motif_search and other algorithms:
#   - score_consensus : sum of hamming distances from the consensus column
#   - score_entropy   : sum of Shannon entropies across columns
# Both take a list of motifs (list[str]) and return a float (lower = better).
# ---------------------------------------------------------------------------

def score_consensus(motifs: list) -> int:
    """
    Scores a set of motifs as the total number of mismatches with the
    column-wise consensus (sum of mistmatches from the max per column).

    Input  : motifs - list of t strings of equal length k
    Output : int - total score (lower = more conserved)

    Complexity: O(t · k)
    """
    t = len(motifs)
    k = len(motifs[0])
    total = 0
    for j in range(k):
        counts = np.zeros(4, dtype=int)
        for i in range(t):
            counts[MAPPING[motifs[i][j]]] += 1
        total += t - counts.max()
    return int(total)


def score_entropy(motifs: list) -> float:
    """
    Scores a set of motifs as the sum of Shannon entropies across columns.
    H(col) = -Σ p_i · log2(p_i)  (0 · log2(0) defined as 0).

    Input  : motifs - list of t strings of equal length k
    Output : float - total entropy (lower = more conserved)

    Complexity: O(t · k)
    """
    t = len(motifs)
    k = len(motifs[0])
    total = 0.0
    for j in range(k):
        counts = np.zeros(4)
        for i in range(t):
            counts[MAPPING[motifs[i][j]]] += 1
        probs = counts / t
        col_entropy = -sum(p * np.log2(p) for p in probs if p > 0)
        total += col_entropy
    return total


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def hamming_distance(s1: str, s2: str) -> int:
    """
    Computes the Hamming distance between two strings of equal length.

    Input  : s1, s2 - strings of equal length
    Output : int - number of positions where characters differ

    Complexity: O(k)
    """
    return sum(a != b for a, b in zip(s1, s2))


def profile_matrix(motifs: list) -> np.ndarray:
    """
    Builds a (4 × k) profile matrix with Laplace pseudocounts
    (each cell initialized to 1 before counting).

    Input  : motifs - list of t strings of equal length k
    Output : np.ndarray shape (4, k), dtype float - probability matrix

    Complexity: O(t · k)
    """
    k = len(motifs[0])
    t = len(motifs)
    p_matrix = np.ones((4, k), dtype=float)
    for motif in motifs:
        for j, base in enumerate(motif):
            p_matrix[MAPPING[base]][j] += 1
    return p_matrix / (t + 4)


def profile_most_probable_kmer(dna_seq: str, k: int, p_matrix: np.ndarray) -> str:
    """
    Returns the k-mer in dna_seq with the highest probability under p_matrix.

    Input  : dna_seq  - DNA string of length n
             k        - k-mer length
             p_matrix - (4 × k) profile matrix
    Output : str - most probable k-mer of length k

    Complexity: O(n · k)
    """
    max_prob = -1.0
    best_kmer = ''
    for i in range(len(dna_seq) - k + 1):
        kmer = dna_seq[i: i + k]
        prob = 1.0
        for j, base in enumerate(kmer):
            prob *= p_matrix[MAPPING[base], j]
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer

def consensus_from_motifs(motifs: list) -> str:
    """
    Returns the motif consensus based on a set of same lenght motifs.

    Input  : motifs - list of k-lenght strings
    
    Output : str - consensus k-mer of lenght k based on most-probable bases

    Complexity: O(t · k)
    """
    k = len(motifs[0])
    consensus = ""
    for i in range(k):
        counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        for motif in motifs:
            counts[motif[i]] += 1
        consensus += max(counts, key=counts.get)
    return consensus

# ---------------------------------------------------------------------------
# BA2b — Median String
# ---------------------------------------------------------------------------

def number_to_pattern(index: int, k: int) -> str:
    """
    Converts an integer index (0 to 4^k - 1) to its corresponding DNA k-mer
    using the lexicographic ordering A=0, C=1, G=2, T=3.

    Input  : index - int in [0, 4^k)
             k     - k-mer length
    Output : str - DNA k-mer of length k

    Complexity: O(k)
    """
    bases = ['A', 'C', 'G', 'T']
    pattern = ""
    for _ in range(k):
        pattern = bases[index % 4] + pattern
        index //= 4
    return pattern


def median_string(k: int, dna_sequences: list) -> str:
    """
    Finds the k-mer that minimizes the total Hamming distance to the set of
    DNA sequences (i.e., the sum of minimum per-sequence distances).

    Exhaustively evaluates all 4^k possible k-mers.

    Input  : k             - int, k-mer length
             dna_sequences - list of t DNA strings
    Output : str - median k-mer of length k

    Complexity: O(4^k · n · k · t)
    """
    min_count = float('inf')
    median = ''
    
    for i in range(4**k):
        pattern_i = number_to_pattern(i, k)
        curr_count = 0
        
        for seq in dna_sequences:
            min_hd = float('inf')
            for j in range(len(seq) - k + 1):
                hd = hamming_distance(pattern_i, seq[j: j+k])
                if hd < min_hd:
                    min_hd = hd
                    if min_hd == 0:
                        break
            curr_count += min_hd
            if curr_count >= min_count:
                break
        
        if curr_count < min_count:
            min_count = curr_count
            median = pattern_i
    
    return median


# ---------------------------------------------------------------------------
# BA2e — Greedy Motif Search with Pseudocounts
# ---------------------------------------------------------------------------

def greedy_motif_search(dna:list, k: int, scoring=score_consensus) -> list:
    """
    Greedy motif search with Laplace pseudocounts. Iterates over all k-mers
    in the first sequence as seeds, then greedily extends the motif set by
    choosing the profile-most-probable k-mer from each subsequent sequence.

    Does not guarantee an optimal solution for the Subtle Motif Problem.

    Input  : k       - int, k-mer length
             dna     - list of t DNA strings of equal length n
             scoring - callable(list[str]) -> float, scoring function to
                       minimize. Options: score_consensus (default),
                       score_entropy.
                       
    Output : tuple(list[str], int) - (best motifs found, their score)

    Complexity: O(n · t · k · n)  →  O(n² · t · k)
    """
    t = len(dna)
    best_motifs = [seq[:k] for seq in dna]

    for w in range(len(dna[0]) - k + 1):
        motifs = [dna[0][w: w + k]]
        for i in range(1, t):
            profile = profile_matrix(motifs)
            motifs.append(profile_most_probable_kmer(dna[i], k, profile))
        if scoring(motifs) < scoring(best_motifs):
            best_motifs = motifs

    return best_motifs, scoring(best_motifs)


# ---------------------------------------------------------------------------
# BA2f — Randomized Motif Search
# ---------------------------------------------------------------------------

def _motifs_from_profile(p_matrix: np.ndarray, dna: list) -> list:
    """
    Given a profile matrix, returns the most probable k-mer from each sequence.

    Input  : p_matrix - (4 × k) profile matrix
             dna      - list of t DNA strings
    Output : list[str] - one k-mer per sequence

    Complexity: O(t · n · k)
    """
    k = p_matrix.shape[1]
    return [profile_most_probable_kmer(seq, k, p_matrix) for seq in dna]


def randomized_motif_search(dna: list, k: int, scoring=score_consensus) -> tuple:
    """
    Single run of the randomized motif search (Monte Carlo).
    Starts from a random set of k-mers (one per sequence) and iteratively
    improves by rebuilding the profile and re-selecting motifs until no
    improvement is found.

    Input  : dna - list of t DNA strings of equal length n
             k   - int, k-mer length
    Output : tuple(list[str], int) - (best motifs found, their score) on this run

    Complexity: O(n · t · k) per iteration; number of iterations varies.
    """
    n = len(dna[0])
    t = len(dna)
    motif_selection = [dna[i][r:r + k]
                       for i, r in ((i, random.randint(0, n - k)) for i in range(t))]

    best_motifs = list(motif_selection)

    while True:
        profile = profile_matrix(motif_selection)
        motif_selection = _motifs_from_profile(profile, dna)
        if scoring(motif_selection) < scoring(best_motifs):
            best_motifs = list(motif_selection)
        else:
            return best_motifs, scoring(best_motifs)


# ---------------------------------------------------------------------------
# BA2g — Gibbs Sampler
# ---------------------------------------------------------------------------

def _profile_random_kmer(dna_i: str, k: int, p_matrix: np.ndarray) -> str:
    """
    Samples a k-mer from dna_i according to the probability distribution
    induced by p_matrix (rather than taking the argmax).

    Input  : dna_i    - DNA string
             k        - k-mer length
             p_matrix - (4 × k) profile matrix
    Output : str - sampled k-mer

    Complexity: O(n · k)
    """
    kmers = [dna_i[j: j + k] for j in range(len(dna_i) - k + 1)]
    probs = np.array([
        float(np.prod([p_matrix[MAPPING[base]][j] for j, base in enumerate(km)]))
        for km in kmers
    ])
    probs /= probs.sum()
    return kmers[np.random.choice(len(kmers), p=probs)]


def gibbs_sampler(dna: list, k: int, N: int, scoring=score_consensus) -> tuple:
    """
    Gibbs sampler for motif finding. At each of N iterations, randomly
    removes one motif, rebuilds the profile from the remaining t-1 motifs,
    and probabilistically samples a new motif for the excluded sequence.

    Input  : dna - list of t DNA strings of equal length n
             k   - int, k-mer length
             N   - int, number of iterations
    Output : tuple(list[str], int) — (best motifs found, their score)

    Complexity: O(N · t · n · k)
    """
    t = len(dna)
    motifs = [seq[r:r + k]
              for seq, r in ((seq, random.randint(0, len(dna[0]) - k)) for seq in dna)]

    best_motifs = list(motifs)

    for _ in range(N):
        i = random.randint(0, t - 1)
        reduced = [motifs[j] for j in range(t) if j != i]
        p_matrix = profile_matrix(reduced)
        motifs[i] = _profile_random_kmer(dna[i], k, p_matrix)
        if scoring(motifs) < scoring(best_motifs):
            best_motifs = list(motifs)

    return best_motifs, scoring(best_motifs)
