#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1k
Title: Generate the Frequency Array of a String
URL: https://rosalind.info/problems/BA1k/

Description:
This script receives a string (DNA sequence) and an integer k.
It returns the frequency of the exact match of all the posible k-mers.
The element position i in the array represents the i k-mer
lexicographically organized, with i = 0 being the first k-mer (AAA...A).

Author: Santiago Wilders Azara 
Date: 2026
"""

import sys

def get_val(symbol):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return mapping.get(symbol.upper(), 0)

def patternToNumber(pattern):
    number = 0
    for char in pattern:
        number = (number * 4) + get_val(char)
    return number

def numberToPattern(index, k):
    bases = ['A', 'C', 'G', 'T']
    pattern = ""
    for _ in range(k):
        pattern = bases[index % 4] + pattern
        index //= 4
    return pattern

def compute_frequencies_rolling(text, k):
    """
    Calcula frecuencias de k-mers en un texto usando Rolling Hash.
    """
    frecuencias = [0] * (4**k)
    n = len(text)
    if n < k: return frecuencias

    # Primer k-mer calculado de forma completa
    current_idx = patternToNumber(text[0:k])
    frecuencias[current_idx] += 1
    
    # Precalculamos el valor para remover el primer nucleótido (4^(k-1))
    mask = 4**(k - 1)

    # El resto se calcula "reciclando" el índice anterior
    for i in range(1, n - k + 1):
        # Rolling Hash: Quitar primero, desplazar, añadir nuevo
        current_idx = ((current_idx % mask) * 4) + get_val(text[i + k - 1])
        frecuencias[current_idx] += 1

    return frecuencias

def main():
    if len(sys.argv) != 2:
        print("Use: python ComputeFrequencies.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            dna_seq = lines[0].strip()
            k= lines[1].strip()
            
            d_match_set = compute_frequencies_rolling(dna_seq, int(k))
            
            print(" ".join(map(str, d_match_set)))
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()