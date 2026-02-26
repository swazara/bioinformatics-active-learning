#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1j
Title: Frequent Words with Mismatches and Reverse Complements Problem
URL: https://rosalind.info/problems/ba1j/

Description:
This script receives a string (DNA sequence) as well as integers k and d.
It then finds and returns the most frequent pattern of length k, taking into
account both the complementary versions of the pattern and with up to d mismatches.

Author: Santiago Wilders Azara 
Date: 2026
"""

import sys

TRANS = str.maketrans("ACGT", "TGCA")
def reverse_complement(dna):
    return dna.translate(TRANS)[::-1]


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
    

def frequent_pattern_mismatches(dna_seq, k, d, n=1):
    counts = {}
    # 1. Conteo de k-meros exactos
    for i in range(len(dna_seq) - k + 1):
        pattern = dna_seq[i : i+k]
        counts[pattern] = counts.get(pattern, 0) + 1
    
    # 2. Construcción del mapa de frecuencias (incluyendo RC)
    frequency_map = {}
    for pattern, count in counts.items():
        neighbors = neighborhood(pattern, d)
        for neighbor in neighbors:
            nc = reverse_complement(neighbor)
            # Sumamos el conteo del patrón original a sus vecinos y sus RC
            frequency_map[neighbor] = frequency_map.get(neighbor, 0) + count
            frequency_map[nc] = frequency_map.get(nc, 0) + count

    # 3. Ordenar por frecuencia y devolver los top n
    # Convertimos el diccionario en una lista de tuplas (patrón, frecuencia)
    # y la ordenamos de mayor a menor frecuencia.
    sorted_patterns = sorted(frequency_map.items(), key=lambda item: item[1], reverse=True)
    
    return sorted_patterns[:n]

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