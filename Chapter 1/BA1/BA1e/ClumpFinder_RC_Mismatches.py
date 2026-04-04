#!/usr/bin/env python3

import sys
from collections import defaultdict
from collections import deque

def reverse_complement(pattern):
    tabla = str.maketrans('ACGT', 'TGCA')
    return pattern.translate(tabla)[::-1]

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

def clumpFinder_with_mismatches_and_rc(secuencia, k, L, t, d=0, use_rc=False):
    n = len(secuencia)
    if n < L: 
        return []

    freq_map = {}
    clumps = {} # Cambiado a dict para guardar la tupla solo la primera vez que alcanza 't'
    memo_neighborhood = {}
    
    def get_neighborhood(kmer):
        if kmer not in memo_neighborhood:
            memo_neighborhood[kmer] = neighborhood(kmer, d)
        return memo_neighborhood[kmer]

    def get_all_variants(kmer):
        variants = set(get_neighborhood(kmer))
        if use_rc:
            rc_kmer = reverse_complement(kmer)
            variants.update(get_neighborhood(rc_kmer))
        return variants

    # 1. Inicialización de la primera ventana de tamaño L
    for i in range(L - k + 1):
        kmer = secuencia[i : i + k]
        for variant in get_all_variants(kmer):
            if variant not in freq_map:
                # 'q' guarda tuplas (indice_en_secuencia, es_exacto)
                freq_map[variant] = {'q': deque(), 'exact': 0, 'var': 0}
            
            is_exact = (variant == kmer)
            freq_map[variant]['q'].append((i, is_exact))
            
            if is_exact:
                freq_map[variant]['exact'] += 1
            else:
                freq_map[variant]['var'] += 1

            if (freq_map[variant]['exact'] + freq_map[variant]['var']) >= t:
                if variant not in clumps:
                    first_pos = freq_map[variant]['q'][0][0]
                    clumps[variant] = (variant, first_pos, freq_map[variant]['exact'], freq_map[variant]['var'])

    # 2. Movimiento de la ventana L por toda la secuencia
    for i in range(1, n - L + 1):
        kmer_out = secuencia[i - 1 : i - 1 + k]
        kmer_in  = secuencia[i + L - k : i + L]
        idx_in   = i + L - k

        # Si entra y sale lo mismo, nada cambia (el 'if' fue removido intencionalmente 
        # para forzar la actualización de los índices en la cola deque)
        
        # --- K-MER SALIENTE ---
        for variant in get_all_variants(kmer_out):
            _, is_exact = freq_map[variant]['q'].popleft()
            if is_exact:
                freq_map[variant]['exact'] -= 1
            else:
                freq_map[variant]['var'] -= 1
                
        # --- K-MER ENTRANTE ---
        for variant in get_all_variants(kmer_in):
            if variant not in freq_map:
                freq_map[variant] = {'q': deque(), 'exact': 0, 'var': 0}
                
            is_exact = (variant == kmer_in)
            freq_map[variant]['q'].append((idx_in, is_exact))
            
            if is_exact:
                freq_map[variant]['exact'] += 1
            else:
                freq_map[variant]['var'] += 1
            
            # Verificamos si se convirtió en un clump
            if (freq_map[variant]['exact'] + freq_map[variant]['var']) >= t:
                if variant not in clumps:
                    first_pos = freq_map[variant]['q'][0][0]
                    clumps[variant] = (variant, first_pos, freq_map[variant]['exact'], freq_map[variant]['var'])

    return list(clumps.values())

def main():
    if len(sys.argv) != 2:
        print("Usage: python ClumpFinderRC&Mismatches.py <input_file.txt>")
        return

    file_path = sys.argv[1]

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()    
            sequence = lines[0].strip()
            k, l_window, t = map(int, lines[1].split())
            
            result_clumps = clumpFinder_with_mismatches_and_rc(sequence, k, l_window, t, d = 1, use_rc=False)

            for tupla in result_clumps:
                patron, pos, exactos, variantes = tupla
                print(f"{patron} (Primera pos: {pos} | Exactos: {exactos} | Variantes/RC: {variantes})")
            
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()