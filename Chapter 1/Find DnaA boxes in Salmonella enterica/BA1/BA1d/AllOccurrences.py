#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1D
Title: Find All Occurrences of a Pattern in a String
URL: https://rosalind.info/problems/ba1d/

Description:
This script finds all starting positions where a given pattern appears as a 
substring within a larger DNA sequence. It accounts for overlapping occurrences.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys

def find_pattern_occurrences(sequence, pattern):
    """
    Finds all starting indices of 'pattern' within 'sequence'.
    Includes overlapping matches.
    """
    indices = []
    pattern_len = len(pattern)

    for i in range(len(sequence) - pattern_len + 1):
        if sequence[i : i + pattern_len] == pattern:
            indices.append(i)
    
    return indices

def main():
    if len(sys.argv) != 2:
        print("Usage: python AllOccurrences.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            # BA1D input format: 
            # Line 1: Pattern
            # Line 2: Genome sequence
            lines = [line.strip() for line in f.readlines() if line.strip()]

            pattern = lines[0]
            sequence = lines[1]
            
            occurrences = find_pattern_occurrences(sequence, pattern)

            # Expected output
            print(" ".join(map(str, occurrences)))
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()