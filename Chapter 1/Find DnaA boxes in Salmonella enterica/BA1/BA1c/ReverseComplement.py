#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1C
Title: Find the Reverse Complement of a String
URL: https://rosalind.info/problems/ba1c/

Description:
This script returns the reverse complement of a DNA sequence.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys

def get_reverse_complement(pattern): 
    """Returns the reverse complement of a given DNA pattern."""
    complementary_nuc = {"A": "T", "T": "A", "C": "G", "G": "C"}
    complement = []
    
    for n in pattern:
        complement.append(complementary_nuc[n.upper()])
    
    reverse_complement = "".join(reversed(complement))
    
    return reverse_complement


def main():
    # Command line argument check
    if len(sys.argv) != 2:
        print("Usage: python ReverseComplement.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            sequence = f.readlines()[0].strip()

            print(get_reverse_complement(sequence))
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' does not exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()