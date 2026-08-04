#!/usr/bin/env python3
"""
Rosalind Problem ID: BA3b
Title: Reconstruct a String from its Genome Path
URL: https://rosalind.info/problems/BA3b/

Description:
Reconstructs a string based on the k-1 equal values between each immediate lines.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys


def string_spelled_from_kmers(kmers):
    if not kmers:
        return
    
    n = len(kmers) 
    k = len(kmers[0])
    string_spelled = kmers[0]

    for i in range(1, n):
        string_spelled += kmers[i][-1]

    return string_spelled
    


def main():
    if len(sys.argv) < 2:
        print("Usage: python StringSpelled.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            kmers = [line.strip() for line in f.readlines() if line.strip()]

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    if kmers:
        print(string_spelled_from_kmers(kmers))

if __name__ == "__main__":
    main()