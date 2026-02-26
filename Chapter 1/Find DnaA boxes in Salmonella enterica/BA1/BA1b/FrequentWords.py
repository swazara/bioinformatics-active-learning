#!/usr/bin/env python3
"""
Rosalind Problem ID: BA1B
Title: Find the Most Frequent Words in a String
URL: https://rosalind.info/problems/ba1b/

Description:
This script identifies the most frequent k-mers (patterns of length k) in a 
given DNA sequence. It builds a frequency table and returns all patterns 
that achieve the maximum count.

Author: Santiago Wilders Azara
Date: 2026
"""

import sys

def get_frequency_table(text, k):
    """Generates a dictionary mapping each k-mer to its occurrence count."""
    frequency_map = {}
    for i in range(len(text) - k + 1):
        pattern = text[i : i+k]
        frequency_map[pattern] = frequency_map.get(pattern, 0) + 1
    
    return frequency_map

def find_frequent_words(text, k):
    """Returns a set of the most frequent k-mers in the text."""
    frequent_patterns = set() 
    frequency_map = get_frequency_table(text, k)
    
    if not frequency_map:
        return frequent_patterns

    # Find the maximum frequency value
    max_freq = max(frequency_map.values())
    
    # Collect all patterns that appear max_freq times
    for pattern, count in frequency_map.items():
        if count == max_freq:
            frequent_patterns.add(pattern)
    
    return frequent_patterns


def main():
    if len(sys.argv) < 2:
        print("Usage: python FrequentWords.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            # BA1B input format: Line 1 = Sequence, Line 2 = k
            sequence = lines[0].strip()
            k = int(lines[1].strip())
            
            print(f"--- Processing: {file_path} (k={k}) ---")
            
            if len(sequence) < k:
                print("Error: The sequence is shorter than k.")
                return
                
            result_patterns = find_frequent_words(sequence, k)
            print(" ".join(result_patterns))
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()