#ifndef GRAPHTRAVERSAL_H
#define GRAPHTRAVERSAL_H
#include "BloomUtil.hpp"
#include "GraphUtil.hpp"
#include "FringeLine.hpp"
#include <set>
#include <fstream>
#include "KmerUtil.hpp"
#include "Utility.hpp"
//static const int min_contig_length=(2*kmer_size+1); //can not use extern variable with constant
//static const uint64_t max_contig_length=10000000;
static const char bin2nt[4]={'A','C','G','T'};
class traverse_kmers
{
    private:
    bloom_util &bloomB;
    branching_kmers *branching_obj;
    uint64_t max_contig_length;
    int max_depth;
    int max_breadth;
    static const int consensus_smilarity=90;
    public:
    traverse_kmers(bloom_util &given_bloomB,branching_kmers *given_branching_obj,
                   uint64_t given_max_len,int given_depth, int given_breadth):bloomB(given_bloomB),
                   branching_obj(given_branching_obj),max_contig_length(given_max_len),
                   max_depth(given_depth),max_breadth(given_breadth){}
    bool find_starting_kmer(kmercode_length branching_kmer, kmercode_length &starting_kmer);
    bool get_path_extension_seed(kmercode_length branching_kmer, kmercode_length &starting_kmer);
    bool explore_branching(kmercode_length start_kmer,int start_strand,char* consensus,
                           int &cons_length,kmercode_length previous_kmer);
    bool explore_branching(kmercode_length start_kmer,int start_strand,char* consensus,
                           int &cons_length,kmercode_length previous_kmer,std::set<kmercode_length> *involved_extensions);
    int find_end_of_branching(kmercode_length start_kmer,int start_strand,kmercode_length &end_kmer,int &end_strand,
                           kmercode_length previous_kmer,std::set<kmercode_length> *involved_extensions);
    std::set<std::string> all_consensus_between(kmercode_length start_kmer,int start_strand,kmercode_length end_kmer,
                                                int end_strand,int traversal_depth,bool &success);
    std::set<std::string> all_consensus_between(kmercode_length start_kmer,int start_strand,kmercode_length end_kmer,
                                                int end_strand,int traversal_depth,std::set<kmercode_length> visited_kmers,
                                                std::string current_consensus,bool &success);
    bool consensus_validation(std::set<std::string> consensus_sequences,char* result, int &result_length);
    bool consensuses_almost_similar(std::set<std::string> consensus_sequences);
    void mark_extensions(std::set<kmercode_length> *involved_extensions);
    int  move_step_forward(kmercode_length current_kmer,int current_strand,bool first_extension,
                           char* nt_new,kmercode_length previous_kmer);
    int  move_step_forward_simple_path(kmercode_length current_kmer,int current_strand,bool first_extension,char* nt_new);
    int traverse(kmercode_length start_kmer,std::vector<char> &contig_sequence,int current_strand,kmercode_length previous_kmer=0);
    int extensions(kmercode_length kmer,int strand,int &nt);
    std::vector< std::pair<int,int> > bubbles_positions;//record start and end positions of traversed bubbles from the latest traverse call()
    void init(uint64_t &assembly_size,uint64_t& number_contigs,uint64_t& max_length);

};
#endif // GRAPHTRAVERSAL_H
