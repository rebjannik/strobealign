#ifndef STROBEALIGN_REFS_HPP
#define STROBEALIGN_REFS_HPP

#include <cstdint>
#include <string>
#include <stdexcept>
#include <numeric>
#include <vector>
#include "exceptions.hpp"
static unsigned char seq_nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

class References {
    typedef std::vector<unsigned int> ref_lengths;
    typedef std::vector<std:: string> ref_names;

public:
    References() { }
    References(
        std::vector<std::string>& string_sequences,
        ref_names names_
    ) : names(std::move(names_)) {

        if (sequences.size() != names.size()) {
            throw std::invalid_argument("lengths do not match");
        }

        //ToDo
        // Here is the conversion function for the byte size
        // Each vector is a byte (or 4 characters), 
        lengths.reserve(sequences.size());

        for (auto& seq : sequences) {
            lengths.push_back(seq.size());
        }

        _total_length = std::accumulate(this->lengths.begin(), this->lengths.end(), (size_t)0);
        
        for ( int i = 0; i < _total_length; ++i){
            for (char nucleotide : sequences[i]) {
                int c = seq_nt4_table[(uint8_t)nucleotide];
                sequences[i] = c < 4 ? nucleotide : 3;
            }
        }
    }

    void add(std::string&& name, std::vector<int>&& sequence);

    static References from_fasta(const std::string& filename);

    size_t size() const {
        return sequences.size();
    }

    size_t total_length() const {
        return _total_length;
    }


    std::vector<std :: vector<int>> sequences;
    ref_names names;
    ref_lengths lengths;
private:
    size_t _total_length;
};

void to_uppercase(std::string& s);

#endif
