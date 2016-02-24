#ifndef GRAPHUTIL_H
#define GRAPHUTIL_H
#include "BloomUtil.hpp"
#include "KmerUtil.hpp"
#include "BinaryStore.hpp"
#include <map>
#include <queue>
#include <set>
extern int min_abun;
typedef std::map<kmercode_length, unsigned int> hash_table;
class branching_kmers
{
private:
    binary_store &solid_kmers;

public:
    hash_table &ht;
    bloom_util &bloomB;
    branching_kmers(binary_store &given_solid_kmers,bloom_util &given_bloomB,hash_table &given_ht)
        :solid_kmers(given_solid_kmers),bloomB(given_bloomB),ht(given_ht) {}
    void compute_branching_kmers(uint64_t &nb_branching_kmers);
    void compute_branching_kmers(kmercode_length kmercode,uint64_t &nb_branching_kmers);
    //void filter_branching_kmers();
    bool is_branching(kmercode_length kmercode);
    unsigned char branching_record(kmercode_length kmercode);
    std::map<kmercode_length, unsigned int>::iterator itr;
    void start_iterator();
    bool next_iterator();
    bool next_branching_node(kmercode_length &kmercode);
    void clear_kmer_abundnce();
    bool is_branching_node(kmercode_length kmercode);
    void mark_node(kmercode_length kmercode);
    void mark_node(kmercode_length given_kmercode, int nt, int strand);
    bool is_marked_branching(kmercode_length kmercode);
    bool is_marked(kmercode_length kmercode);
    bool is_marked(kmercode_length kmercode,char nt,int strand);
    bool is_empty_();
    int get_batch_of_kmers(kmercode_length *kmers_codes,int max_batch_size);


};

#endif // GRAPHUTIL_H
