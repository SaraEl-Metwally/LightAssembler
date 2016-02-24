#include "GraphUtil.hpp"
int min_abun;
unsigned char branching_kmers::branching_record(kmercode_length kmercode)
{
    unsigned char record=0;
    kmercode_length new_kmercode=0;
    int nt,strand;
    for(nt=0; nt<4; nt++)
    {
        //forward extension --------> strand 0
        strand=0;
        new_kmercode=next_kmer(kmercode,nt,&strand);
        if(bloomB.contain(new_kmercode))
            record|=1<<nt;
        // <------------ backward extension strand 1
        strand=1;
        new_kmercode=next_kmer(kmercode,nt,&strand);
        if(bloomB.contain(new_kmercode))
            record|=1<<nt+4;

    }
    return record;
}
bool branching_kmers::is_branching(kmercode_length kmercode)
{
    int num_fw_links=0;
    int num_bw_links=0;
    int i;
    unsigned char record=branching_record(kmercode);
    for(i=0; i<4; i++)
        num_fw_links += (record>>i)&1;
    for(i=4; i<8; i++)
        num_bw_links += (record>>i)&1;

    return !(num_fw_links==1 && num_bw_links==1);
}
/*
void branching_kmers::filter_branching_kmers()
{
    //min_abun=3;
    std::cout<< "--- minimum abundance: "<< min_abun <<std::endl;
    uint64_t nb_branching_kmers=0;
    itr = ht.begin();
    while (itr != ht.end())
    {
        if(itr->second<= min_abun)
        {
            ht.erase(itr++);
        }
        else
        {
            itr->second=0;//No Need for abundance any more.
            ++itr;
            nb_branching_kmers++;
        }
    }
    //you can store it in vector now because the value int does not need to store
    // because all kmers now have the default value of min_abun
    std::cout<<"--- Number of filtered branching kmers: "<<nb_branching_kmers<<std::endl;

}
*/
void branching_kmers::start_iterator()
{
    itr = ht.begin();

}
bool branching_kmers::next_iterator()
{
    if(itr==ht.end()) return false;
    return true;

}
bool branching_kmers::next_branching_node(kmercode_length &given_kmercode)
{
    //you should verify the result of this method.
    if(next_iterator())
    {
        given_kmercode= itr->first ;
        itr++;
        return true;
    }
    return false;
}


int branching_kmers::get_batch_of_kmers(kmercode_length *kmers_codes,int max_batch_size)
{

   int temp = solid_kmers.read_element_buffer(kmers_codes,max_batch_size);
   
   return temp;

}

void branching_kmers::compute_branching_kmers(kmercode_length kmercode,uint64_t &nb_branching_kmers)
{

       if(is_branching(kmercode))
        {
            if (ht.find(kmercode) == ht.end() )
            {
                ht[kmercode]=1;
                nb_branching_kmers++; //number of unique nb_branching kmers without repetition.
            }

            else
                ht[kmercode]=ht[kmercode]+1;
        }


}

void branching_kmers::compute_branching_kmers(uint64_t &nb_branching_kmers)
{
    solid_kmers.rewind_all();
    kmercode_length kmercode=0;

    //Note: the file of solid kmers not sorted, so it contains many and many repeated kmers.


    while(solid_kmers.read_element(&kmercode))
    {
        if(is_branching(kmercode))
        {
            if (ht.find(kmercode) == ht.end() )
            {
                ht[kmercode]=1;
                nb_branching_kmers++; //number of unique nb_branching kmers without repetition.
            }

            else
                ht[kmercode]=ht[kmercode]+1;

        }
        
    }

}
void branching_kmers::clear_kmer_abundnce()
{
    itr = ht.begin();
    while (itr != ht.end())
    {
        itr->second=0;
        itr++;
    }
}

bool branching_kmers::is_empty_()
{
    if(ht.size()==0)
        return true;
    else
        return false;

}

bool branching_kmers::is_branching_node(kmercode_length given_kmercode)
{
    if (ht.find(given_kmercode) == ht.end())
      return false;
      else
      return true;

}
void branching_kmers::mark_node(kmercode_length given_kmercode)
{

    if(is_branching_node(given_kmercode))//if it is indexed in a ht
    {
        unsigned int val=ht[given_kmercode];// previous value contains status of 8 possible extension
        val|= 1<<8; //9 bits 0-->7 8 bits for 8 possible extension + 1 bit for mark as used in the assembly graph or not.
        ht[given_kmercode]= val;
        // you can use flag to assert marking as minia
        //flag=true;

    }
    for(int strand=0;strand<2;strand++)
    {
        for (int nt=0;nt<4;nt++)
        {
            int neighbor_str=strand;
            // Note: the following method returns the minimum of the two neighbor codes (code and revcomp) and
            // set the nighbor_str with the value of chosen strand.
            kmercode_length neighbor_kmr=next_kmer(given_kmercode,nt,&neighbor_str);

            if(!bloomB.contain(neighbor_kmr))
              continue;
            if(!is_branching_node(neighbor_kmr))//mark only branching nodes or their neighbors of branching ones, no simple nodes
              continue;
            int neighbor_rev_str=1-neighbor_str;
            int nt_code;
            if(strand==0)
                nt_code=rev_nt2int(code2nt(given_kmercode,0));
            else
                nt_code=code2nt(given_kmercode,kmer_size-1);

            mark_node(neighbor_kmr,nt_code,neighbor_rev_str);
        }

    }

   //return flag;

}
//record branching kmer info..
void branching_kmers::mark_node(kmercode_length given_kmercode, int nt, int strand)
{
    if(!is_branching_node(given_kmercode)) //if it is not branching node, ignore it. //if it is indexed in a ht
        return;// false;
    if(is_marked(given_kmercode,nt,strand))
    return;
    unsigned int val=0;
    val=ht[given_kmercode];
    if(strand==0)
        val|=1<<(nt);
    else
        val|=1<<(nt+4);
    ht[given_kmercode]=val;
}

bool branching_kmers::is_marked_branching(kmercode_length given_kmercode)
{
    unsigned int val=ht[given_kmercode];
    return (val&(1<<8)) !=0; //test if it is marked

}
bool branching_kmers::is_marked(kmercode_length given_kmercode)
{
    if(is_branching_node(given_kmercode)) //if it is indexed in a ht
        return is_marked_branching(given_kmercode);
    for(int strand=0;strand<2;strand++)
    {
        for (int nt=0;nt<4;nt++)
        {
            int neighbor_str=strand;
            // Note: the following method returns the minimum of the two neighbor codes (code and revcomp) and
            // set the nighbor_str with the value of chosen strand.
            kmercode_length neighbor_kmr=next_kmer(given_kmercode,nt,&neighbor_str);
            if(!bloomB.contain(neighbor_kmr))
              continue;
            if(!is_branching_node(neighbor_kmr))//mark only branching nodes or their neighbors of branching ones, no simple nodes
              continue;
            int neighbor_rev_str=1-neighbor_str;
            int nt_code;
            if(strand==0)
                nt_code=rev_nt2int(code2nt(given_kmercode,0));
            else
                nt_code=code2nt(given_kmercode,kmer_size-1);

            if(is_marked(neighbor_kmr,nt_code,neighbor_rev_str))
                return true;
        }

    }

  return false;
}
bool branching_kmers::is_marked(kmercode_length given_kmercode,char nt,int strand)
{
    if(!is_branching_node(given_kmercode))//kmer not indexed
        return false;
    unsigned int val=ht[given_kmercode];
    int extend_nt_mrk;
    if(strand==0)
        extend_nt_mrk=(val>>nt)&1;
    else
        extend_nt_mrk=(val>>(nt+4))&1;
  return extend_nt_mrk==1;
}
