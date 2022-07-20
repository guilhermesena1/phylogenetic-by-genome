#include <iostream>
#include <string>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <stdio.h>

#define VERBOSE

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::ifstream;
using std::copy;
using std::ostream_iterator;
using std::back_inserter;
using std::endl;
using std::runtime_error;
static const uint8_t encode_char[256] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //4
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //17
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //33
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //49
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //@,A-O
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //P-Z
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, //`,a-o
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //p-z
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};


struct KmerStats {
  KmerStats() { kmer_count = vector<uint32_t>(num_kmers, 0); }
  void count_kmers(const string &chrom);
  vector<uint32_t> kmer_count;
  static const uint32_t kmer_size = 12;
  static const uint32_t num_kmers = (1 << (2*kmer_size));
};

inline void
shift_hash_key(const uint8_t c, uint32_t &mer) {
  static const uint32_t hash_mask = KmerStats::num_kmers - 1;
  mer = ((mer << 2) | encode_char[c]) & hash_mask;
}

void
KmerStats::count_kmers(const string &chrom) {
  uint32_t mer = 0;

  auto itr(begin(chrom));
  const auto lim(end(chrom));
  for (size_t i = 0; itr != lim && i < kmer_size - 1; ++i, ++itr)
    shift_hash_key(*itr, mer);

  for (; itr != lim; ++itr) {
    shift_hash_key(*itr, mer);
    ++kmer_count[mer];
  }
}

inline void
process_chrom(string &chrom, KmerStats &v) {
  // makes sure letters are only ACGT:
  for (auto it(begin(chrom)); it != end(chrom); ++it)
    *it = ((encode_char[static_cast<uint8_t>(*it)] == 4) ?  "ACGT"[rand()%4] : *it);

  v.count_kmers(chrom);
  chrom.clear();
}

void
process_species(const string &file, KmerStats &v) {
  // get file size
  ifstream in(file);

  if (!in)
    throw runtime_error("cannot open file " + file);

  static const size_t RESERVE_SIZE = 250000000;
  string chrom;
  chrom.reserve(RESERVE_SIZE);

  string line;
  while (getline(in, line)) {
    if (line[0] != '>')
      copy(begin(line), end(line), back_inserter(chrom));
    else {
      process_chrom(chrom, v);
    }
  }

  process_chrom(chrom, v); // last chromosome after EOF
  in.close();
}

int main(int argc, const char **argv) {
  if (argc != 2) {
    cout << "usage: ./phylo <input-species.txt>" << endl;
    return 0;
  }

  // ensures the k-mer size used fits in a 32-bit number
  static_assert(KmerStats::kmer_size <= 16);

  vector<string> species;
  vector<string> files;

  string tmp1, tmp2;
  ifstream in(argv[1]);
  while (in >> tmp1 >> tmp2) {
    species.push_back(tmp1);
    files.push_back(tmp2);
  }
  in.close();
  omp_set_num_threads(8);

  copy(begin(species), end(species), std::ostream_iterator<string>(cout, "\t"));
  cout << "\n";
  // print the headers: the species names, separated by tabs
  vector<KmerStats> v(species.size());
#pragma omp parallel for
  for (size_t i = 0; i < species.size(); ++i) {
#ifdef VERBOSE

#pragma omp critical
    {
      cerr << "processing " << species[i] << "...\n";
    }
#endif
    process_species(files[i], v[i]);
  }

#ifdef VERBOSE
  cerr << "writing output\n";
#endif
  const size_t num_species = species.size();
  for (size_t i = 0; i < KmerStats::num_kmers; ++i) {
    for (size_t j = 0; j < num_species; ++j)
      printf("%d\t", v[j].kmer_count[i]);
    printf("\n");
  }
  return 0;
}
