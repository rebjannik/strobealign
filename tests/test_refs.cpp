#include <fstream>
#include "zstr.hpp"
#include "doctest.h"
#include "refs.hpp"

TEST_CASE("References::add") {
    References references;
    references.add(std::string("thename"), std::string("ACGT"));

    CHECK(references.names.size() == 1);
    CHECK(references.names[0] == "thename");
    //CHECK(references.sequences[0] == "ACGT");
    CHECK(references.substr(0, 0, 4) == "ACGT");
    CHECK(references.lengths[0] == 4);
}

TEST_CASE("References::from_fasta") {
    auto references = References::from_fasta("tests/phix.fasta");
    CHECK(references.names.size() == 1);
    CHECK(references.names[0] == "NC_001422.1");
    CHECK(references.sequences.size() == 1);
    CHECK(references.lengths.size() == 1);
    CHECK(references.lengths[0] == 5386);
    CHECK(references.total_length() == 5386);
}

TEST_CASE("References::from_fasta parse error") {
    REQUIRE_THROWS_AS(References::from_fasta("tests/phix.1.fastq"), InvalidFasta);
}

TEST_CASE("Reference FASTA not found") {
    REQUIRE_THROWS_AS(References::from_fasta("does-not-exist.fasta"), InvalidFasta);
}

TEST_CASE("Reference uppercase") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs
            << ">ref1\n"
            << "acgt\n\n"
            << ">ref2\n"
            << "aacc\ngg\n\ntt\n"
            << ">empty\n"
            << ">empty_at_end_of_file";
    }
    auto refs = References::from_fasta("tmpref.fasta");
    std::remove("tmpref.fasta");
    CHECK(refs.sequences.size() == 2);
    //CHECK(refs.sequences[0].size() == 4);
    CHECK(refs.sequences[0].size() == 1);
    // CHECK(refs.sequences[0] == "ACGT");
    CHECK(refs.substr(0,0,4) == "ACGT");
    //CHECK(refs.sequences[1].size() == 8);
    CHECK(refs.sequences[1].size() == 2);
    // CHECK(refs.sequences[1] == "AACCGGTT");
    CHECK(refs.substr(1, 0, 8) == "AACCGGTT");
    CHECK(refs.names.size() == 2);
    CHECK(refs.lengths.size() == 2);
}

TEST_CASE("Reference gzipped") {
    {
        zstr::ofstream ofs("tmpref.fasta.gz");
        ofs << ">ref1\n"
            << "acgt\n";
    }
    auto refs = References::from_fasta("tmpref.fasta.gz");
    std::remove("tmpref.fasta.gz");
    CHECK(refs.sequences.size() == 1);
    CHECK(refs.names[0] == "ref1");
    //CHECK(refs.sequences[0] == "ACGT");
    CHECK(refs.substr(0, 0, 4) == "ACGT");
}

TEST_CASE("Invalid contig name") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs << ">name\b\n";
    }
    REQUIRE_THROWS_AS(References::from_fasta("tmpref.fasta"), InvalidFasta);
    std::remove("tmpref.fasta");
}

TEST_CASE("Invalid contig name") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs << ">\0x80\nACGT\n";
    }
    REQUIRE_THROWS_AS(References::from_fasta("tmpref.fasta"), InvalidFasta);
    std::remove("tmpref.fasta");
}

TEST_CASE("Tab and space in FASTA header") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs
            << ">ab c\nACGT\n"
            << ">de\tf\nACGT\n";
    }
    auto refs = References::from_fasta("tmpref.fasta");
    std::remove("tmpref.fasta");
    CHECK(refs.size() == 2);
    CHECK(refs.names[0] == "ab");
    CHECK(refs.names[1] == "de");
}

TEST_CASE("Duplicate contig name") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs
            << ">abc\nAAAA\n"
            << ">abc\nCCCC\n";
    }
    REQUIRE_THROWS_AS(References::from_fasta("tmpref.fasta"), InvalidFasta);
    std::remove("tmpref.fasta");
}

TEST_CASE("Accept some special characters in contig name") {
    {
        std::ofstream ofs("tmpref.fasta");
        ofs
            << "><>;abc\nAAAA\n"
            << ">abc\nCCCC\n";
    }
    auto refs = References::from_fasta("tmpref.fasta");

    CHECK(refs.size() == 2);
    CHECK(refs.names[0] == "<>;abc");
    std::remove("tmpref.fasta");
}
