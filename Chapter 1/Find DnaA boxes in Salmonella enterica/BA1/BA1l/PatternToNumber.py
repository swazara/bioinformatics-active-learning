#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1l
Title: Implement PatternToNumber
URL: https://rosalind.info/problems/BA1l/

Description:
This script receives a dna pattern and returns its corresponding number
using a 4-base with each number realtive to its lexicographical order.

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

"""""
def numberToPattern(index, k):
    bases = ['A', 'C', 'G', 'T']
    pattern = ""
    for _ in range(k):
        pattern = bases[index % 4] + pattern
        index //= 4
    return pattern
"""

def main():
    if len(sys.argv) != 2:
        print("Use: python PatternToNumber.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            pattern = lines[0].strip()
            
            print(patternToNumber(pattern))
            
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()