#ifndef BloomUtil_H
#define BloomUtil_H
#include "BloomFilter.hpp"
#include <stdint.h>
class bloom_util
{
private:
    uint64_t bf_size;
    bloom_parameters bf_para;
    bloom_filter bf;
    int nb_threads;
public:
    bloom_util(double fprate=0.01): bf_size(10000003),bf_para(10000003,fprate), bf(bf_para),nb_threads(1){}
    bloom_util(uint64_t bf_size,double fprate=0.01): bf_size(bf_size),bf_para(bf_size,fprate), bf(bf_para),nb_threads(1){}
    double get_occupancy()
    {
        return bf.occupancy();
    }

    double get_false_positive_rate()
    {
        return bf.GetActualFP();
    }
    
    #ifdef largeintlib
    int add(kmercode_length &val,bool flag=false)
    {
        if(nb_threads>=1&&flag&&(bf.contains(val)))
          return 0;
        bf.insert(val);
         return 1;
    }
    #endif
    #ifdef _LP64
    int add(__uint128_t &val,bool flag=false)
    {
        if(nb_threads>=1&&flag&&(bf.contains(val)))
          return 0;
        bf.insert(val);
          return 1;
    }
    #endif
    int add(uint64_t &val,bool flag=false)
    {
        if(nb_threads>=1&&flag&&(bf.contains(val)))
          return 0;
        bf.insert(val);
          return 1;
    }
    #ifdef largeintlib
    bool contain(LargeInt<kmer_precision> &val)
    {
        return bf.contains(val)	;
    }
    #endif
    #ifdef _LP64
    bool contain(__uint128_t &val)
    {
        return bf.contains(val) ;

    }
    #endif
    bool contain(uint64_t &val)
    {
        return bf.contains(val) ;

    }

   void set_nb_threads( int nb_t) 
    { 
		nb_threads = nb_t ;
		bf.SetNumOfThreads( nb_t ) ;
    }

};
#endif
