# Author: Santiago Wilders Azara
# Date: 2026
# Description: Bioinformatics utility functions for k-mer analysis and DNA sequence processing.
# Developed with assistance from Claude.

"""
bio_utils.py — Bioinformatics utility functions.

Notation used in complexity comments
--------------------------------------
  n  : length of the DNA sequence
  k  : length of a k-mer / pattern
  L  : window size
  t  : minimum occurrence threshold
  d  : maximum number of mismatches allowed
  N  : size of the d-neighborhood of a k-mer  →  N = Σ_{i=0}^{d} C(k,i)·3^i
"""

from collections import defaultdict, deque
# Optional dependencies — only required for visualization functions
try:
    import seaborn as sns
    import matplotlib.pyplot as plt
except ImportError:
    sns = None
    plt = None

# ---------------------------------------------------------------------------
# Core primitives
# ---------------------------------------------------------------------------
def hamming_distance(s1: str, s2: str) -> int:
    """
    Computes the Hamming distance between two strings of equal length.
    
    Complexity: O(k)
    """
    return sum(a != b for a, b in zip(s1, s2))

def reverse_complement(pattern: str) -> str:
    """
    Returns the reverse complement of a DNA string.

    Complexity: O(k)
    """
    tabla = str.maketrans("ACGT", "TGCA")
    return pattern.translate(tabla)[::-1]


def neighborhood(pattern: str, d: int) -> set:
    """
    Returns the set of all DNA strings within Hamming distance d of *pattern*.

    Strategy: recursive generation — at each position either keep the original
    nucleotide or (if mismatches remain) substitute it with any of the 3 others.

    Complexity: O(k · N)  where N = Σ_{i=0}^{d} C(k,i)·3^i
    """
    def generate(idx: int, remaining_d: int):
        if idx == len(pattern):
            return [""]
        results = []
        original = pattern[idx]
        for suffix in generate(idx + 1, remaining_d):
            results.append(original + suffix)
        if remaining_d > 0:
            for base in "ACGT":
                if base != original:
                    for suffix in generate(idx + 1, remaining_d - 1):
                        results.append(base + suffix)
        return results

    return set(generate(0, d))


# ---------------------------------------------------------------------------
# G-C Skew functions
# ---------------------------------------------------------------------------

def skew_list(dna_seq: str, step: int = 1) -> list:
    """
    Computes the cumulative G-C skew array for a DNA sequence.
    Position 0 is always 0; position i reflects the skew after base i*step.

    step : sample every `step` nucleotides (default 1 = full resolution).
           Useful for visualization of large genomes. Does not affect
           skew_min(), which always runs at full resolution.

    Complexity: O(n)
    """
    skew_values = {"G": 1, "C": -1, "A": 0, "T": 0}
    curr = 0
    history = [0]
    for i, base in enumerate(dna_seq.upper()):
        curr += skew_values.get(base, 0)
        if (i + 1) % step == 0:
            history.append(curr)
    return history
    
def skew_min(dna_seq: str) -> list:
    """
    Returns all positions (1-based) where the skew reaches its global minimum.
    Position 0 is included in the candidates since the sequence may start there.

    Complexity: O(n)
    """
    skew_values = {"G": 1, "C": -1, "A": 0, "T": 0}
    curr, min_val = 0, 0
    positions = [0]
    for i, base in enumerate(dna_seq.upper(), 1):
        curr += skew_values.get(base, 0)
        if curr < min_val:
            min_val = curr
            positions = [i]
        elif curr == min_val:
            positions.append(i)
    return positions


# ---------------------------------------------------------------------------
# Pattern counting
# ---------------------------------------------------------------------------

def count_pattern_occurrences(sequence: str, pattern: str) -> tuple:
    """
    Counts exact occurrences of *pattern* in *sequence* (with overlaps).

    Returns: (count, [starting_positions])

    Complexity: O(n · k)
    """
    k = len(pattern)
    count = 0
    positions = []
    for i in range(len(sequence) - k + 1):
        if sequence[i: i + k] == pattern:
            count += 1
            positions.append(i)
    return count, positions


# ---------------------------------------------------------------------------
# Frequent patterns with mismatches (grouped by canonical / RC pair)
# ---------------------------------------------------------------------------

def frequent_patterns_with_mismatches(dna_seq: str, k: int, d: int, n: int = 10) -> list:
    """
    Finds the top-*n* most frequent k-mers in *dna_seq* allowing up to *d*
    mismatches.  Reverse-complement pairs are grouped together and reported
    under their lexicographically smaller representative.

    Returns a list of tuples sorted by total score (descending):
        (canonical_kmer, total_score, score_forward, score_reverse_complement)

    Complexity: O(n · N + 4^k · k)
      Dominant term is usually O(n · N) for short k with small d.
    """
    # Step 1 — exact counts
    # O(n)
    exact_counts: dict = {}
    for i in range(len(dna_seq) - k + 1):
        kmer = dna_seq[i: i + k]
        exact_counts[kmer] = exact_counts.get(kmer, 0) + 1

    # Step 2 — propagate to neighbors (forward only; RC grouping happens in step 3)
    # O(n · N)
    mismatch_scores: dict = {}
    for kmer, count in exact_counts.items():
        for neighbor in neighborhood(kmer, d):
            mismatch_scores[neighbor] = mismatch_scores.get(neighbor, 0) + count

    # Step 3 — group by RC pair using canonical (lexicographically smaller) form
    # O(4^k · k)
    grouped: list = []
    seen: set = set()

    for pattern in sorted(mismatch_scores):
        if pattern in seen:
            continue
        rc = reverse_complement(pattern)
        score_fwd = mismatch_scores.get(pattern, 0)
        score_rev = mismatch_scores.get(rc, 0)

        if pattern == rc:                  # palindromic k-mer
            total = score_fwd
        else:
            total = score_fwd + score_rev

        # Canonical = lexicographically smaller of the pair
        canonical = min(pattern, rc)
        if canonical == pattern:
            grouped.append((canonical, total, score_fwd, score_rev))
        else:
            # forward strand is the RC, so swap the directional scores
            grouped.append((canonical, total, score_rev, score_fwd))

        seen.add(pattern)
        seen.add(rc)

    grouped.sort(key=lambda x: x[1], reverse=True)
    return grouped[:n]


# ---------------------------------------------------------------------------
# Clump finding  (generalised: supports mismatches and reverse complement)
# ---------------------------------------------------------------------------

def find_clumps(
    sequence: str,
    k: int,
    L: int,
    t: int,
    d: int = 0,
    use_rc: bool = False,
) -> list:
    """
    Finds all k-mers that form an (k, L, t)-clump in *sequence*.

    A k-mer forms a clump when it appears at least *t* times (counting
    d-mismatches and, optionally, its reverse complement) within some
    window of length *L*.

    Returns a list of tuples:
        (kmer, first_window_position, exact_count, variant_count)

    Algorithm: sliding window of width L.
    - On each slide, remove the leftmost k-mer's contribution and add the
      newly entering k-mer's contribution.
    - Neighborhoods are memoised to avoid recomputation.

    Complexity: O(n · N)  where N is the neighborhood size.
    """
    n = len(sequence)
    if n < L:
        return []

    memo_neighborhood: dict = {}

    def get_variants(kmer: str) -> set:
        if kmer not in memo_neighborhood:
            variants = neighborhood(kmer, d)
            if use_rc:
                variants.update(neighborhood(reverse_complement(kmer), d))
            memo_neighborhood[kmer] = variants
        return memo_neighborhood[kmer]

    # freq_map[variant] = {'q': deque of (position, is_exact), 'exact': int, 'var': int}
    freq_map: dict = {}
    clumps: dict = {}   # variant -> first tuple that reached threshold

    def add_kmer(kmer: str, pos: int):
        for variant in get_variants(kmer):
            if variant not in freq_map:
                freq_map[variant] = {"q": deque(), "exact": 0, "var": 0}
            is_exact = variant == kmer
            freq_map[variant]["q"].append((pos, is_exact))
            if is_exact:
                freq_map[variant]["exact"] += 1
            else:
                freq_map[variant]["var"] += 1
            if (freq_map[variant]["exact"] + freq_map[variant]["var"]) >= t:
                if variant not in clumps:
                    first_pos = freq_map[variant]["q"][0][0]
                    clumps[variant] = (
                        variant, first_pos,
                        freq_map[variant]["exact"],
                        freq_map[variant]["var"],
                    )

    def remove_kmer(kmer: str):
        for variant in get_variants(kmer):
            _, is_exact = freq_map[variant]["q"].popleft()
            if is_exact:
                freq_map[variant]["exact"] -= 1
            else:
                freq_map[variant]["var"] -= 1

    # Initialise first window
    for i in range(L - k + 1):
        add_kmer(sequence[i: i + k], i)

    # Slide window
    for i in range(1, n - L + 1):
        remove_kmer(sequence[i - 1: i - 1 + k])
        idx_in = i + L - k
        add_kmer(sequence[idx_in: idx_in + k], idx_in)

    return list(clumps.values())


def clean_clump_noise(clump_results: list, d: int = 1) -> list:
    """
    From a list of clump tuples produced by find_clumps, groups k-mers
    that are within Hamming distance d of each other and keeps only the
    best representative (highest exact count, then highest total count).

    Complexity: O(c² · k)  where c = len(clump_results)
    """
    remaining = list(clump_results)
    result = []

    while remaining:
        # Take the current best candidate (highest exact count)
        remaining.sort(key=lambda x: (x[2], x[2] + x[3]), reverse=True)
        best = remaining.pop(0)
        result.append(best)

        # Remove all k-mers within Hamming distance d of the best
        remaining = [
            entry for entry in remaining
            if hamming_distance(entry[0], best[0]) > d
        ]

    return sorted(result, key=lambda x: x[1])

# ---------------------------------------------------------------------------
# Results visualization
# ---------------------------------------------------------------------------
    
def highlight_sequence_html(
    sequence: str,
    patterns: list,
    colors: list = None,
    offset: int = 0,
    line_width: int = 60,
    markers: list = None,
    font_size: int = 13,
) -> None:
    """
    Renders a DNA sequence as colored HTML in a Jupyter notebook.

    Parameters
    ----------
    sequence   : DNA string to display
    patterns   : list of patterns to highlight. Each entry can be:
                   - a string              → exact match only (d=0)
                   - a tuple (str, int)    → pattern with d mismatches allowed
    colors     : list of hex color strings, one per pattern.
                 If None, a default palette is used.
                 If fewer colors than patterns, cycles through the list.
    offset     : genome coordinate of sequence[0], used for position labels
    line_width : characters per line (default 60)
    markers    : list of (absolute_position, label, color) tuples to mark
                 specific landmarks (e.g. skew minimum) with a vertical tick
    font_size  : display font size in px (default 13)
    """
    from IPython.display import display, HTML

    DEFAULT_COLORS = [
        "#86efac",  # green
        "#fdba74",  # orange
        "#93c5fd",  # blue
        "#f9a8d4",  # pink
        "#c084fc",  # purple
        "#fde68a",  # yellow
    ]

    # 1. Normalize patterns → list of (pattern, d)
    normalized = []
    for entry in patterns:
        if isinstance(entry, tuple):
            normalized.append(entry)
        else:
            normalized.append((entry, 0))

    # 2. Assign colors
    palette = colors if colors is not None else DEFAULT_COLORS
    pattern_colors = [palette[i % len(palette)] for i in range(len(normalized))]

    # 3. Build position → color map (last pattern in list wins on overlap)
    n = len(sequence)
    pos_colors = {}

    for (pattern, d), color in zip(normalized, pattern_colors):
        neighbors = neighborhood(pattern, d)
        plen = len(pattern)
        for i in range(n - plen + 1):
            if sequence[i: i + plen] in neighbors:
                for j in range(plen):
                    pos_colors[i + j] = color

    # 4. Convert markers to relative positions
    marker_positions = {}
    if markers:
        for abs_pos, label, color in markers:
            rel = abs_pos - offset
            if 0 <= rel < n:
                marker_positions[rel] = (label, color)

    # 5. Build HTML
    html = (
        f"<pre style='font-family: monospace; font-size: {font_size}px; "
        f"line-height: 2; white-space: pre-wrap; word-break: break-all;'>"
    )

    for i, char in enumerate(sequence):
        if i % line_width == 0:
            html += (
                f"\n<span style='color: gray; font-size: {font_size - 2}px;'>"
                f"{i + offset:>8,}  </span>"
            )
        if i in marker_positions:
            label, mcolor = marker_positions[i]
            html += (
                f"<span style='color: {mcolor}; font-weight: bold;' "
                f"title='{label}'>▼</span>"
            )
        if i in pos_colors:
            html += f"<span style='background-color: {pos_colors[i]};'>{char}</span>"
        else:
            html += char

    html += "</pre>"

    # 6. Legend
    legend = (
        f"<div style='font-family: monospace; margin-bottom: 6px; "
        f"font-size: {font_size}px;'>"
    )
    for (pattern, d), color in zip(normalized, pattern_colors):
        label = f"{pattern} (d≤{d})" if d > 0 else pattern
        legend += (
            f"<span style='background-color: {color}; padding: 2px 8px; "
            f"margin-right: 8px;'>{label}</span>"
        )
    if markers:
        for _, label, mcolor in markers:
            legend += (
                f"<span style='border-left: 3px solid {mcolor}; padding-left: 6px; "
                f"margin-right: 8px;'>{label}</span>"
            )
    legend += "</div>"

    display(HTML(legend + html))

def plot_gc_skew(
    sequence: str,
    step: int = 500,
    title: str = "G-C Skew",
    figsize: tuple = (14, 5),
    color: str = "darkcyan",
    marker_color: str = "crimson",
) -> tuple:
    """
    Plots the G-C skew along a DNA sequence and marks the minimum position(s).

    Parameters
    ----------
    sequence     : DNA string
    step         : sampling resolution for visualization (default 500)
    title        : plot title
    figsize      : figure size tuple
    color        : skew line color
    marker_color : color for the minimum position marker

    Returns
    -------
    min_positions : list of positions (1-based) where skew is minimal
    min_value     : the minimum skew value reached
    """
    
    if sns is None or plt is None:
        raise ImportError("seaborn and matplotlib are required for visualization. Install them with: pip install seaborn matplotlib")
    
    # 1. Compute skew
    skew_vals      = skew_list(sequence, step=step)
    real_positions = [i * step for i in range(len(skew_vals))]

    # 2. Find minimum
    min_positions = skew_min(sequence)
    min_value     = min(skew_vals)
    min_pos       = min(min_positions)

    # 3. Plot
    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=figsize)

    sns.lineplot(x=real_positions, y=skew_vals, color=color, linewidth=1.2, ax=ax, label="G-C Skew")

    ax.axvline(x=min_pos, color=marker_color, linewidth=1.2, linestyle="--", alpha=0.8)
    ax.annotate(
        f"min = {min_value}\npos = {min_pos:,}",
        xy=(min_pos, min_value),
        xytext=(min_pos + len(sequence) * 0.02, min_value + 5),
        fontsize=9,
        color=marker_color,
        arrowprops=dict(arrowstyle="->", color=marker_color, lw=1.0),
    )

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Position in genome (bp)", fontsize=12)
    ax.set_ylabel("Skew value (G − C)", fontsize=12)
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.show()

    print(f"Minimum skew value  : {min_value}")
    print(f"Minimum position(s) : {[f'{p:,}' for p in min_positions]}")

    return min_positions, min_value