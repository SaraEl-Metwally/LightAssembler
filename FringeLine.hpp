#ifndef FRINGELINE_H
#define FRINGELINE_H
#include "BloomUtil.hpp"
#include "GraphUtil.hpp"
#include "KmerUtil.hpp"
#include <queue>
#include <set>
struct node
{
    kmercode_length kmer;
    int strand;
    int nt;
    node(kmercode_length given_kmer,int given_strand,int given_nt):kmer(given_kmer),strand(given_strand),nt(given_nt){}
    //check
    bool operator<(const node &other) const
    {
        if(kmer !=other.kmer)
            return (kmer< other.kmer);
        if(strand !=other.strand)
            return (strand<other.strand);
        return (nt<other.nt);

    }
};
class fringe_line
{   //check names
 private:
    kmercode_length start_kmer;//initialize
    kmercode_length previous_kmer;
    int start_strand;//initialize
    bloom_util &bloomB;
    branching_kmers *branching_obj;
    std::set<kmercode_length> *involved_extensions;
    std::set<kmercode_length> already_queued_nodes;
    std::queue<node> fringe;
public:
    bool in_branching;
    int depth;
    fringe_line(kmercode_length start_kmer,
               int start_strand,
               bloom_util &bloomB,
               branching_kmers *branching_obj,
               std::set<kmercode_length> *involved_extensions,
               kmercode_length previous_kmer=0,
               bool in_branching=true);
    bool next_depth();
    int current_breadth();
    node peak();
    bool check_in_branching(kmercode_length kmer, int strand);
};
#endif // FRINGELINE_H
