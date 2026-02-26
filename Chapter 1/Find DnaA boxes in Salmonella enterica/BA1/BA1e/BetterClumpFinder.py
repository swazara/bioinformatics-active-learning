# En teoría un mejor algoritmo para encontrar clumps, el problema, Python.

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

def ComputeFrequenciesRolling(text, k):
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

def clumpFinder(secuencia, valores):
    k, L, t = map(int, valores.split())
    n = len(secuencia)
    if n < L: return []

    clump_marks = [0] * (4**k)
    mask = 4**(k - 1)

    # 1. Inicialización del primer bloque L
    frecuencias = ComputeFrequenciesRolling(secuencia[0:L], k)
    
    for i in range(4**k):
        if frecuencias[i] >= t:
            clump_marks[i] = 1

    # 2. Inicializamos los dos punteros de Rolling Hash
    # index_start: el k-mer que empieza en la posición i
    # index_end: el k-mer que empieza en la posición (i + L - k)
    index_start = patternToNumber(secuencia[0:k])
    index_end = patternToNumber(secuencia[L-k : L])

    # 3. Movimiento de la ventana L por toda la secuencia (O(n))
    for i in range(1, n - L + 1):
        # --- K-MER SALIENTE ---
        frecuencias[index_start] -= 1
        
        # Prepara el index_start para la próxima iteracón
        index_start = ((index_start % mask) * 4) + get_val(secuencia[i + k - 1])
        
        # --- K-MER ENTRANTE ---
        index_end = ((index_end % mask) * 4) + get_val(secuencia[i + L - 1])
        
        frecuencias[index_end] += 1
        
        # Verificamos si este nuevo k-mer es un clump
        if frecuencias[index_end] >= t:
            clump_marks[index_end] = 1

    # 4. Conversión final
    return [numberToPattern(i, k) for i in range(4**k) if clump_marks[i] == 1]

def main():
    if len(sys.argv) != 2:
        print("Uso: python script.py archivo.txt")
        return

    archivo = sys.argv[1]
    try:
        with open(archivo, 'r') as f:
            lineas = f.read().splitlines()
            if len(lineas) < 2: return
            secuencia = lineas[0].strip()
            valores = lineas[1].strip()
            
            resultado = clumpFinder(secuencia, valores)
            print(" ".join(resultado))
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()