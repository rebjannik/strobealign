#ifndef STROBEALIGN_REFS_HPP
#define STROBEALIGN_REFS_HPP

#include <cstdint>
#include <string>
#include <stdexcept>
#include <numeric>
#include <vector>
#include "exceptions.hpp"

extern unsigned char seq_nt4_table[256]; 

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
            std::vector<uint8_t> packed_seq;
            pack_sequence(str, packed_seq);
            sequences.push_back(std::move(packed_seq));
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
        return unpack_sequence(sequences[ref_id], lengths[ref_id]);
    }

    std::string unpack_sequence(const std::vector<uint8_t>& seq, size_t length) const {
        std::string result;
        result.reserve(length);
        size_t idx = 0;

        for (uint8_t byte : seq) {
                for (int shift = 6; shift >= 0 && idx < length; shift -= 2, idx++) {
                    uint8_t val = (byte >> shift) & 0x3;
                    result += "ACGT"[val];
                }
            }
        return result;
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
            uint8_t byte = sequences[ref_id][i / 4];
            int shift = 6 - 2 * (i % 4);
            uint8_t val = (byte >> shift) & 0x3;
            result += "ACGT"[val];
        }
        return result;
    }

    std::vector<std::vector<uint8_t>> sequences;
    ref_names names;
    ref_lengths lengths;

    private:
    size_t _total_length;

    void pack_sequence(const std::string& str, std::vector<uint8_t>& seq) const {
        uint8_t byte = 0;
        int count = 0;

        for (char c : str){
            uint8_t val = seq_nt4_table[static_cast<unsigned char>(c)];
            byte = (byte << 2) | val;
            count++;
            
            if (count == 4){
                seq.push_back(byte);
                byte = 0;
                count = 0;
            }
        }

        if (count > 0) {
            byte <<= (2 * (4 - count));
            seq.push_back(byte);
        }
    }

};



void to_uppercase(std::string& s);

#endif
