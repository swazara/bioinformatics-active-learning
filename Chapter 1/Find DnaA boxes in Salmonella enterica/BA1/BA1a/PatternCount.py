#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1A
Title: Compute the Number of Times a Pattern Appears in a Text
URL: https://rosalind.info/problems/ba1a/

Description:
Basic pattern-matching algorithm to count the total occurrences of a 
specific k-mer within a larger DNA sequence. 

Author: Santiago Wilders Azara
Date: 2026
"""

import sys

def count_pattern_occurrences(sequence, pattern):
    """Counts the number of times a pattern appears in a sequence, including overlaps."""
    count = 0
    pattern_len = len(pattern)
    for i in range(len(sequence) - pattern_len + 1):
        if sequence[i : i + pattern_len] == pattern:
            count += 1
    
    return count


def main():
    if len(sys.argv) < 2:
        print("Usage: python PatternCount.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            # BA1A input format: Line 1 = Sequence, Line 2 = Pattern
            sequence = lines[0].strip()
            pattern = lines[1].strip()
            
            print(f"--- Processing: {file_path} (Pattern = {pattern}) ---")
            
            total_count = count_pattern_occurrences(sequence, pattern)
            print(f"{total_count} occurrences found.")
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()