#include "GraphTraversal.hpp"
void traverse_kmers::mark_extensions(std::set<kmercode_length> *involved_extensions)
{
    if(branching_obj==NULL)
        return;
    for(std::set<kmercode_length>::iterator it=involved_extensions->begin();it!=involved_extensions->end();++it)
    {
        branching_obj->mark_node(*it);
    }


}
bool traverse_kmers::get_path_extension_seed(kmercode_length branching_kmer, kmercode_length &starting_kmer)
{
    for(int strand =0;strand <2;strand++)
    {
     for(int nt=0;nt<4;nt++)
     {
         int current_str=strand;
         kmercode_length current_kmer=next_kmer(branching_kmer,nt,&current_str);
         if(bloomB.contain(current_kmer))
         {
             if(branching_obj->is_branching_node(current_kmer))
             continue;
             if(branching_obj->is_marked(current_kmer)) //marked before
             continue;
             branching_obj->mark_node(current_kmer);
             starting_kmer=current_kmer;
             return true;
         }

     }

    }
   return false;
}
int traverse_kmers::find_end_of_branching(kmercode_length start_kmer,int start_strand,kmercode_length &end_kmer,int &end_strand,
                                          kmercode_length previous_kmer,std::set<kmercode_length> *involved_extensions)
{
    bool in_branching = true;
    fringe_line fringline(start_kmer,start_strand,bloomB,branching_obj,involved_extensions,previous_kmer,in_branching);
    do
    {
        bool go_on =fringline.next_depth();
        if(! go_on)
            return 0;
        if(fringline.depth> max_depth)
            return 0;
        if(fringline.current_breadth()>max_breadth)
            return 0;
        if(fringline.current_breadth()==0)
            return 0;
        if(fringline.current_breadth()==1 &&(branching_obj==NULL || (!branching_obj->is_branching_node(fringline.peak().kmer))))
            break;
    }
    while(1);
    if(fringline.current_breadth()==1)
    {
        node end_node=fringline.peak();
        end_kmer=end_node.kmer;
        end_strand=end_node.strand;
        return fringline.depth;
    }
   return 0;

}
std::set<std::string> traverse_kmers::all_consensus_between(kmercode_length start_kmer,int start_strand,kmercode_length end_kmer,
                                                int end_strand,int traversal_depth,bool &success)
{
  std::set<kmercode_length> visited_kmers;
  visited_kmers.insert(start_kmer);
  std::string current_consensus;
  success=true;
  return all_consensus_between(start_kmer,start_strand,end_kmer,end_strand,
                               traversal_depth,visited_kmers,current_consensus,success);
}
std::set<std::string> traverse_kmers::all_consensus_between(kmercode_length start_kmer,int start_strand,kmercode_length end_kmer,int end_strand,
                                            int traversal_depth,std::set<kmercode_length> visited_kmers,
                                            std::string current_consensus,bool &success)
{

    std::set<std::string> consensus_sequences;
    if(traversal_depth < -1)
    {
        success=false;
        return consensus_sequences;
    }
    if(start_kmer==end_kmer)
    {
      consensus_sequences.insert(current_consensus);
      return consensus_sequences;
    }
    //traverse all neighbors
    for(int nt=0;nt<4;nt++)
    {
        int new_strand=start_strand;
        kmercode_length new_kmer=next_kmer(start_kmer,nt,&new_strand);
        if(bloomB.contain(new_kmer))
        {
            if(visited_kmers.find(new_kmer)!=visited_kmers.end())//Tandem Repeats: Bubbles of Loops.
            {
                success=false;
                return consensus_sequences;
            }
            std::string extended_consensus_seq(current_consensus);
            extended_consensus_seq.append(1,bin2nt[nt]);
            std::set<kmercode_length> used_kmers(visited_kmers);
            used_kmers.insert(new_kmer);
            std::set<std::string> new_consensus_sequences=all_consensus_between(new_kmer,new_strand,end_kmer,end_strand,
                                                          traversal_depth-1,used_kmers,extended_consensus_seq,success);
            consensus_sequences.insert(new_consensus_sequences.begin(),new_consensus_sequences.end());
            if(consensus_sequences.size()> (unsigned int)max_breadth)
                success=false;
        }//end if
        if(success==false) //force stopping because there are too many consensuses reached
            return consensus_sequences;
    }
    return consensus_sequences;

}
bool traverse_kmers::consensuses_almost_similar(std::set<std::string> consensus_sequences)
{
    for(std::set<std::string>::iterator it_a=consensus_sequences.begin();it_a!=consensus_sequences.end();++it_a)
    {
        std::set<std::string>::iterator it_b= it_a;
        advance(it_b,1);
        while(it_b != consensus_sequences.end())
        {
            if(needleman_wunch(*it_a,*it_b)*100 <consensus_smilarity)
               return false;
            advance(it_b,1);
        }

    }

    return true;
}
bool traverse_kmers::consensus_validation(std::set<std::string> consensus_sequences,char* result, int &result_length)
{
    int mean=0;
    int path_number=0;
    //compute mean and stdev of all bubble paths
    for(std::set<std::string>::iterator it=consensus_sequences.begin();it!=consensus_sequences.end();++it)
    {
        mean +=(*it).length();
        path_number++;

    }
    mean /=consensus_sequences.size();
    double stdev=0;
    for(std::set<std::string>::iterator it=consensus_sequences.begin();it!=consensus_sequences.end();++it)
    {
        int consensus_length=(*it).length();
        stdev +=pow(fabs(consensus_length-mean),2);
    }
    stdev =sqrt(stdev/consensus_sequences.size());
    if(mean>max_depth)
        return false; //traverse large bubbles is not allowed.
    if(consensus_sequences.size()==1&&mean >kmer_size+1)//dead_end length should be <k+1
        return false; //traverse large dead_ends is not allowed.
    if(stdev>mean/5)
        return false; // traverse bubbles with paths that don't have roughly the same length is not allowed.
    if(!consensuses_almost_similar(consensus_sequences))
        return false; //check consensus sequences similarity.
    //Now all paths are filtered to choose among them.
    std::string chosen_consensus =*consensus_sequences.begin();
    result_length=chosen_consensus.length();
    if(result_length>max_depth) //chosen consensus is longer than max_depth
        return false;
    chosen_consensus.copy(result,result_length);
    return true;
}
bool traverse_kmers::explore_branching(kmercode_length start_kmer,int start_strand,char* consensus,
                           int &cons_length,kmercode_length previous_kmer)
{
    std::set<kmercode_length> *involved_extensions =new std::set<kmercode_length>;
    bool flag= explore_branching(start_kmer,start_strand,consensus,cons_length,previous_kmer,involved_extensions);
    delete involved_extensions ;
    return flag;

}
//return true if the branching is successfully traversed and marking all involved nodes.
bool traverse_kmers::explore_branching(kmercode_length start_kmer,int start_strand,char* consensus,
                           int &consensus_length,kmercode_length previous_kmer,std::set<kmercode_length> *involved_extensions)
{
  kmercode_length end_kmer=0; //initialized variable
  int end_strand=0; //initialized variable
  int traversal_depth = find_end_of_branching(start_kmer,start_strand,end_kmer,end_strand,previous_kmer,involved_extensions);
  //the previous method will store the end kmer in end_kmer and it is associated end_strand in end_strand.
  if(!traversal_depth)
     return false; // it is a complex bubble.
  std::set<std::string> consensus_sequences;
  bool success=false;//initialized variable
  //find all consensus sequences between start and end nodes.
  consensus_sequences=all_consensus_between(start_kmer,start_strand,end_kmer,end_strand,traversal_depth+1,success);
  if(!success)
    return false;
 //path validation based on sequence similarity.
 bool valid = consensus_validation(consensus_sequences,consensus,consensus_length);
 if(!valid)
    return false;
 //mark traversed nodes.
 mark_extensions(involved_extensions);
 return true;
}
bool traverse_kmers::find_starting_kmer(kmercode_length branching_kmer,kmercode_length &starting_kmer)
{
    int total_depth=0;
    if(!get_path_extension_seed(branching_kmer,starting_kmer))
        return false;
    for(int strand=0;strand<2;strand++)
    {
        kmercode_length previous_kmer=0;
        int previous_strand=0;
        //BFS to verify that this path is not inside a bubble or tip
        fringe_line fringeline(starting_kmer,strand,bloomB,branching_obj,NULL,0,false);//
        do
        {
            bool go_on =fringeline.next_depth();
            if(!go_on)
                break;
            if(fringeline.depth>max_depth || fringeline.current_breadth()>max_breadth)
                break;
            if(fringeline.current_breadth()==0)
                break;
            char consensus[max_depth+1];
            int consensus_length=0;
            if(fringeline.current_breadth()<=1)
            {
                kmercode_length current_kmer=0;
                if(fringeline.current_breadth()==1)
                {
                    node current_node=fringeline.peak();
                    current_kmer=current_node.kmer;
                }

                if((previous_kmer!=0)&& branching_obj->is_branching_node(previous_kmer))
                {
                    std::set<kmercode_length> involved_extensions;
                    branching_kmers *save_branching=branching_obj;
                    branching_obj=NULL;
                    if(explore_branching(previous_kmer,1-previous_strand,consensus,
                                         consensus_length,current_kmer,&involved_extensions))
                    {
                      if(involved_extensions.find(starting_kmer)!=involved_extensions.end())//it is visited before..
                         {
                             branching_obj=save_branching;
                             return false;//starting kmer is in a tip/bubble path starting from current kmer
                         }//if involved
                    }//if explore
                   branching_obj=save_branching;

                }//if previous kmer

            }//if #nodes in fringeline <= 1

            //update previous kmer
            if(fringeline.current_breadth()==1)
            {
                node current_node=fringeline.peak();
                previous_kmer=current_node.kmer;
                previous_strand=current_node.strand;
            }
            else
                previous_kmer=0;

        }//do
        while(1);
        total_depth +=fringeline.depth;

    }//for strand
    if(total_depth<(kmer_size+1))//avoid assemble regions that do not produce long contigs
        return false;
    return true;
}
int traverse_kmers::extensions(kmercode_length kmer,int strand,int &nt)
{
    //examine immediate neighbors only, you can extend this method to examine more deeper paths to detect dead ends.
    int nb_extensions=0;
    for(int n=0;n<4;n++)
    {
      int current_strand =strand;
      kmercode_length current_kmer=next_kmer(kmer,n,&current_strand);
      if(bloomB.contain(current_kmer))
      {
          nt=n;
          nb_extensions++;

      }

    }
    return nb_extensions;


}
int traverse_kmers::move_step_forward_simple_path(kmercode_length current_kmer,int current_strand,bool first_extension,char* nt_new)
{
    int nb_extensions=0;
    int chosen_nt=0;
    nb_extensions=extensions(current_kmer,current_strand,chosen_nt);
    if(nb_extensions==1)
    {
        int second_strand=current_strand;
        kmercode_length second_kmer= next_kmer(current_kmer,chosen_nt,&second_strand);
        int second_nt=0;
        int in_branching_degree=0;
        in_branching_degree =extensions(second_kmer,1-second_strand,second_nt);
        if(in_branching_degree>1)
            return -2;// next_kmer has multiple in-branching paths

        *nt_new=bin2nt[chosen_nt];
        return 1;//good extension is found;


    }
    if(nb_extensions>1) //this kmer has multiple extension paths, it is an out-branching kmer.
        return -1;

    return 0; //if a dead end is reached.

}
int traverse_kmers::move_step_forward(kmercode_length current_kmer,int current_strand,bool
                  first_extension,char* nt_new,kmercode_length previous_kmer)
{
    int simple_path=move_step_forward_simple_path(current_kmer,current_strand,first_extension,nt_new);
    if(simple_path>0)
        return 1; //it is simple path: it is a simple non-branching kmer.
    //bubble exploration..
    int consensus_length=0;
    bool success =explore_branching(current_kmer,current_strand,nt_new,consensus_length,previous_kmer);
    if(!success)
        return 0;
    return consensus_length;

}
int traverse_kmers::traverse(kmercode_length start_kmer,std::vector<char> &contig_sequence,int start_strand,kmercode_length previous_kmer)
{
     kmercode_length current_kmer=start_kmer;
     int current_strand=start_strand;
     int extension_length=0;
     char new_nt[max_depth+1];
     int traverse_status=0;
     bool circular_region=false;
     int bubble_start=0,bubble_end=0;
     bubbles_positions.clear();
     while((traverse_status=move_step_forward(current_kmer,current_strand,extension_length==0,new_nt,previous_kmer)))
     {
            if(traverse_status<0) //-1:out-branching kmer,-2:in-branching kmer, 0:Dead_end
            break;
            if(traverse_status>1)// >1 it is bubble and it represents its length  1: simple path
            bubble_start=extension_length;
            for(int nt=0;nt<traverse_status;nt++)
            {
                contig_sequence[extension_length]=new_nt[nt];
                extension_length++;
                previous_kmer=current_kmer;
                current_kmer=next_kmer(current_kmer,nt2int(new_nt[nt]),&current_strand);
                branching_obj->mark_node(current_kmer);
                if(current_kmer==start_kmer)//circular region
                    circular_region=true;

            }//for inspect bubble...
            if(traverse_status>1)
            {
                bubble_end=extension_length;
                bubbles_positions.push_back(std::make_pair(bubble_start,bubble_end));

            }
            if(circular_region)
                break;
            if(extension_length>max_contig_length)
                break;

     }//end_while
     return extension_length;

}

void traverse_kmers::init(uint64_t &assembly_size,uint64_t& number_contigs,uint64_t& max_length)
{
    kmercode_length branch_code=0;
    uint64_t left_extension_length=0,right_extension_length=0,contig_length=0;
    uint64_t  nb_contigs=0,nb_nts=0,nb_branching_kmers=0,max_contig_len=0,max_left_len=0,max_right_len=0;
    int min_contig_length=(2*kmer_size+1);
    std::vector<char> left_extension(max_contig_length);
    std::vector <char> right_extension(max_contig_length);
    char kmer_chars[kmer_size+1];
    std::ofstream file_assembly(return_file_name(assembly_file).c_str(),std::ios::out);
    if (file_assembly.is_open())
    {
      while(branching_obj->next_branching_node(branch_code))
         {
           kmercode_length start_kmer=0;
           //nb_branching_kmers++;
           //std::cout<<" Branching no: " <<nb_branching_kmers<<std::endl;
           while(find_starting_kmer(branch_code,start_kmer))
             {
                code2seq(start_kmer,kmer_chars);
                right_extension_length=traverse(start_kmer,right_extension,0);
                max_right_len=std::max(right_extension_length,max_right_len);
                left_extension_length=traverse(start_kmer,left_extension,1);
                max_left_len=std::max(left_extension_length,max_left_len);
                revcomp_sequence(left_extension,left_extension_length);
                std::string seq1(left_extension.begin(),left_extension.begin()+left_extension_length);
                std::string start_kmer_chars(kmer_chars,kmer_size);
                std::string seq2(right_extension.begin(),right_extension.begin()+right_extension_length);
                std::string contig_seq=seq1+start_kmer_chars+seq2;
                contig_length=left_extension_length+kmer_size+right_extension_length;
                if(contig_length >= min_contig_length)
                {
                    max_contig_len=std::max(max_contig_len,contig_length);
                    nb_contigs++;
                    file_assembly <<">"<<nb_contigs<<"___length__"<<contig_length<<std::endl;
                    file_assembly <<contig_seq<<std::endl;
                    nb_nts+=contig_length;
                }//Finish Writing*/


              }// kmers on the path*/
              left_extension.clear();
              right_extension.clear();
             // nb_branching_kmers++; //print message of total number of branching kmers.
        }//branching nodes
        file_assembly.close();
        assembly_size=nb_nts;
        number_contigs=nb_contigs;
        max_length=max_contig_len;
        //std::cout<<"--- in this assembly session:   "<<std::endl;
        //std::cout<<"--- total number of nts: "<<nb_nts<<" are assembled into: "<<nb_contigs<<" contigs"<<std::endl;
        //std::cout<<std::endl;
        //std::cout<<"--- maximum left traversal"<<max_left_len<<std::endl;
        //std::cout<<"--- maximum right traversal"<<max_right_len<<std::endl;
    }
    else
    {
        std::cerr << "--- can't open output contig file " <<return_file_name(assembly_file).c_str() << std::endl;
        exit(1);

    }


}//method

