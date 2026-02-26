#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1m
Title: Implement NumberToPattern
URL: https://rosalind.info/problems/BA1m/

Description:
This script receives a number and an integer k, returns its corresponding pattern.


Author: Santiago Wilders Azara 
Date: 2026
"""

import sys


def numberToPattern(index, k):
    bases = ['A', 'C', 'G', 'T']
    pattern = ""
    for _ in range(k):
        pattern = bases[index % 4] + pattern
        index //= 4
    return pattern


def main():
    if len(sys.argv) != 2:
        print("Use: python NumberToPattern.py <input.txt>")
        return

    try:
        with open(sys.argv[1], 'r', encoding='utf-8') as f:
            lines = f.readlines()
            number = lines[0].strip()
            k = lines[1].strip()
            
            print(numberToPattern(int(number), int(k)))
            
            
    except FileNotFoundError:
        print(f"Error: The '{sys.argv[1]}' file doesn't exist.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()