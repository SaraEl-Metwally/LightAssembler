#include "FringeLine.hpp"
//set of nodes having equal depth in the BFS
fringe_line::fringe_line(kmercode_length start_kmer,
                       int start_strand,
                       bloom_util &bloomB,
                       branching_kmers *branching_obj,
                       std::set<kmercode_length> *involved_extensions,
                       kmercode_length previous_kmer,
                       bool in_branching):
                       start_kmer(start_kmer),
                       start_strand(start_strand),
                       bloomB(bloomB),
                       branching_obj(branching_obj),
                       involved_extensions(involved_extensions),
                       previous_kmer(previous_kmer),
                       in_branching(in_branching),
                       depth(0)
 {
      already_queued_nodes.insert(start_kmer);
      already_queued_nodes.insert(previous_kmer);
      node first_node(start_kmer,start_strand,-1);
      fringe.push(first_node);
  }
int fringe_line::current_breadth()
{
  return fringe.size();
}
node  fringe_line::peak()
{
   return fringe.front();
}
bool fringe_line::next_depth()
{
    std::queue<node> new_fringe;
    while(!fringe.empty())
    {
        node current=fringe.front();
        fringe.pop();
        kmercode_length current_kmer =current.kmer;
        int current_strand =current.strand;
        if(in_branching&& check_in_branching(current_kmer,current_strand))//complex bubble: large in-branching inside: longer than 3k
            return false;
        //enqueue all possible extension neighbors except ones that are already in the fringe.
        for(int nt=0;nt<4;nt++)
        {
            kmercode_length new_kmer=current_kmer;
            int new_strand=current_strand;
            int from_nt=(current.nt==-1)? nt:current.nt;
            new_kmer=next_kmer(new_kmer,nt,&new_strand);
            //note &new_strand if it is defined as pointer but passing it by address and it will change there.
            if(already_queued_nodes.find(new_kmer)!=already_queued_nodes.end())
                continue;//this node is already in the fringe
            if(bloomB.contain(new_kmer))
            {
                //is_branching_node: test hash table of branching nodes.
                //is_branching : test the branching structure associated with this node
                //if(branching_obj !=NULL && branching_obj->is_branching(new_kmer))
                if(branching_obj !=NULL && branching_obj->is_branching_node(new_kmer))
                    if(branching_obj->is_marked_branching(new_kmer))
                      return false;
                node new_node(new_kmer,new_strand,from_nt);
                new_fringe.push(new_node);
                already_queued_nodes.insert(new_kmer);
                if(involved_extensions!= NULL)
                involved_extensions->insert(new_kmer);


            }//if

        }//for


    }//while
    fringe=new_fringe;
    ++depth;
    return true;

}
//detect any in-branching longer than 3k
bool fringe_line::check_in_branching(kmercode_length from_kmer, int from_strand)
{
    for(int nt=0;nt<4;nt++)
    {
        int current_strand=1-from_strand;
        kmercode_length current_kmer=next_kmer(from_kmer,nt,&current_strand);
        //exclude kmers that already queued which contains the previously
        //traversed kmers according to the first extension.
        //you should understand this in deep details.
        if(already_queued_nodes.find(current_kmer)!=already_queued_nodes.end())
            continue;
        if(bloomB.contain(current_kmer))
        {
            //check larger in-branching by creating a new FringeLine inside this one.i.e. need to go to deeper depth.
            fringe_line fringeline(current_kmer,current_strand,bloomB,branching_obj,involved_extensions,from_kmer,false);
            do
            {
                bool stop =fringeline.next_depth();
                if(!stop)
                    break;
                if(fringeline.depth>3*kmer_size)//limit depth
                    break;
                if(fringeline.current_breadth()>10)//limit breadth
                    break;
                if(fringeline.current_breadth()==0) //no more in-branching nodes
                    break;
            }
            while(1);

            if(fringeline.current_breadth()>0)
                return true;

        }

    }
    return false;

}
