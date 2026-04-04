
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>

// Complexity: O(k)
int hamming_distance(const std::string& a, const std::string& b) {
    int dist = 0;
    for (size_t i = 0; i < a.size(); ++i)
        dist += (a[i] != b[i]);
    return dist;
}

// Complexity: O(k)
std::string number_to_pattern(int index, int k) {
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string pattern(k, 'A');
    for (int i = k - 1; i >= 0; --i) {
        pattern[i] = bases[index % 4];
        index /= 4;
    }
    return pattern;
}

// Complexity: O(4^k · n · k · t) — with branch-and-bound pruning
std::string median_string(int k, const std::vector<std::string>& dna) {
    int t = dna.size();
    int min_dist = std::numeric_limits<int>::max();
    std::string median = "";
    int total_patterns = 1;
    for (int i = 0; i < k; ++i) total_patterns *= 4;

    for (int i = 0; i < total_patterns; ++i) {
        std::string pattern = number_to_pattern(i, k);
        int curr_dist = 0;

        for (int s = 0; s < t; ++s) {
            const std::string& seq = dna[s];
            int n = seq.size();
            int min_hd = k;

            for (int j = 0; j <= n - k; ++j) {
                int hd = hamming_distance(pattern, seq.substr(j, k));
                if (hd < min_hd) min_hd = hd;
                if (min_hd == 0) break;  // early exit
            }

            curr_dist += min_hd;
            if (curr_dist >= min_dist) goto next_pattern;  // branch-and-bound
        }

        if (curr_dist < min_dist) {
            min_dist = curr_dist;
            median = pattern;
        }

        next_pattern:;
    }
    return median;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./median_string <input_file>\n";
        return 1;
    }

    std::ifstream file(argv[1]);
    if (!file) {
        std::cerr << "Error: cannot open file\n";
        return 1;
    }

    int k;
    file >> k;
    file.ignore();

    std::vector<std::string> dna;
    std::string line;
    while (std::getline(file, line))
        if (!line.empty()) dna.push_back(line);

    std::cout << median_string(k, dna) << "\n";
    return 0;
}