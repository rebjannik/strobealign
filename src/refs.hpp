#ifndef STROBEALIGN_REFS_HPP
#define STROBEALIGN_REFS_HPP

#include <cstdint>
#include <string>
#include <stdexcept>
#include <numeric>
#include <vector>
#include "exceptions.hpp"

extern unsigned char seq_nt4_table[256]; // Declare the array

class References {
    typedef std::vector<unsigned int> ref_lengths;
    typedef std::vector<std:: string> ref_names;

public:
    References() : _total_length(0) { }
    References(
        const std::vector<std::string>& string_sequences,
        const ref_names& names_
    ) : names(names_), _total_length(0) {
        for (const auto& str : string_sequences) {
            std::vector<int> seq;
            seq.reserve(str.size());
            for (char c : str) {
                seq.push_back(seq_nt4_table[static_cast<unsigned char>(c)]);
            }
            sequences.push_back(std::move(seq));
            lengths.push_back(str.length());
            _total_length += str.length();
        }
        if (sequences.size() != names.size()) {
            throw std::invalid_argument("Lengths do not match");
        }
    }
    
    std::string operator[](size_t ref_id) const {
        if (ref_id >= sequences.size()) {
            throw std::out_of_range("Invalid reference ID");
        }
    
        std::string result;
        result.reserve(sequences[ref_id].size());
        for (int i : sequences[ref_id]) {
            result += "ACGT"[i];
        }
        return result;
    }
    
    void convert(const std::string& str, std::vector<int>& seq) const {
        for (char c : str) {
            seq.push_back(seq_nt4_table[static_cast<unsigned char>(c)]);
        }
    }
    
    void add(std::string&& name, std::string&& sequence);

    static References from_fasta(const std::string& filename);

    size_t size() const {
        return sequences.size();
    }

    size_t total_length() const {
        return _total_length;
    }

    std::string substr (size_t ref_id, size_t start, size_t length) const{
        std:: string result;
       
        result.reserve(length);
        for (size_t i = start; i < start + length; ++i) {
            result += "ACGT"[sequences[ref_id][i]];
    }
        return result;
    }

    std::vector<std :: vector<int>> sequences;
    ref_names names;
    ref_lengths lengths;

    private:
    size_t _total_length;
};

void to_uppercase(std::string& s);

#endif
