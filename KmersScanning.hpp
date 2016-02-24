#ifndef KMERSCANNING_H
#define KMERSCANNING_H
#include "BloomUtil.hpp"
#include "KmerUtil.hpp"
#include "MathUtil.hpp"
#include "BinaryStore.hpp"
#include "Utility.hpp"
#include "GraphUtil.hpp"
#include "GraphTraversal.hpp"
#include "ReadsParsing.hpp"
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <typeinfo>


struct threads_arg_parameters_extrapolation
{
  uint64_t total_bases;
  uint64_t total_reads;
  reads_parsing *reads;
  pthread_mutex_t *lock_A ;
  pthread_mutex_t *lock_S ;

};
struct threads_arg_sample_kmers
{
  uint64_t total_sample_kmers;
  uint64_t total_bases;
  uint64_t total_reads;
  bloom_util *sample_kmers;
  reads_parsing *reads;
  pthread_mutex_t *lock_A ;
  pthread_mutex_t *lock_S ;

};
struct threads_arg_store_trusted_kmers
{
  std::vector<int> threshold ;
  uint64_t total_kmers_b;
  uint64_t occured_kmers;
  bloom_util *sample_kmers;
  bloom_util *trusted_kmers;
  binary_store *solid_kmers;
  reads_parsing *reads;
  pthread_mutex_t *lock_F ;
  pthread_mutex_t *lock_S ;
};

struct threads_arg_branching_kmers
{
  uint64_t total_kmers_branching;
  uint64_t total_kmers_file;
  pthread_mutex_t *lock_H ;
  pthread_mutex_t *lock_F ;
  pthread_mutex_t *lock_S ;
  branching_kmers* kmers_branching_obj;
};



void *pass_zero_parameters_extrapolation_multi_threads( void *arg )
{
  struct threads_arg_parameters_extrapolation* pass_zero_arg = (struct threads_arg_parameters_extrapolation*)arg;
  int max_batch_size= READ_BUFFER_PER_THREAD ;
  int curr_batch_size;
  one_read *read_batch= (one_read*)malloc(sizeof(one_read) * max_batch_size) ;
  uint64_t total_bases=0;uint64_t total_reads=0;
  while(1)
  {

         pthread_mutex_lock( pass_zero_arg->lock_A ) ;
         curr_batch_size = pass_zero_arg->reads->get_batch_of_reads(read_batch, max_batch_size);
         pthread_mutex_unlock( pass_zero_arg->lock_A );
         if(curr_batch_size==0)
           break;
         for(int i=0;i<curr_batch_size;i++)
         {
            std::string read_seq(read_batch[i].read_seq);
            total_reads++;
            total_bases+=read_seq.length();
         } //for batch

  }//while
 pthread_mutex_lock( pass_zero_arg->lock_S ) ;
 pass_zero_arg->total_bases=pass_zero_arg->total_bases+total_bases;
 pass_zero_arg->total_reads=pass_zero_arg->total_reads+total_reads;
 pthread_mutex_unlock( pass_zero_arg->lock_S );
 free( read_batch ) ;
 pthread_exit(NULL);
 return NULL ;

}

void pass_one_bloom_a(std::string read,bloom_util *bloomA,int &nb_sample_kmers_p)
{
    int i;
    int j=gap_size;
    kmercode_length kmercode=0,kmercode_can=0;
    for(i=0; i<kmer_size; ++i)
    {

        kmercode=kmercode*4+nt2int(read[i]);


    }
    kmercode_can=get_canonical_kmer_code(kmercode); 
    bloomA->add(kmercode_can);
    nb_sample_kmers_p++;
    for(i=1; i<read.length()-kmer_size+1; ++i)
    {
        kmercode=(kmercode*4+nt2int(read[i+(kmer_size-1)])& kmer_mask);
        
        if(i==j)
         {
           kmercode_can=get_canonical_kmer_code(kmercode);
           bloomA->add(kmercode_can);
           j=j+gap_size;
           nb_sample_kmers_p++;
         }

    }


}//method*/

void *pass_one_bloom_a_multi_threads( void *arg )
{
  struct threads_arg_sample_kmers* pass_one_arg = (struct threads_arg_sample_kmers*)arg;
  int max_batch_size= READ_BUFFER_PER_THREAD ;
  int curr_batch_size; 
  one_read *read_batch= (one_read*)malloc(sizeof(one_read) * max_batch_size) ;
  uint64_t total_kmers_a=0;uint64_t total_bases=0;uint64_t total_reads=0;
  while(1)
  {

       pthread_mutex_lock( pass_one_arg->lock_A ) ;
       curr_batch_size = pass_one_arg->reads->get_batch_of_reads(read_batch, max_batch_size);
       pthread_mutex_unlock( pass_one_arg->lock_A );
       if(curr_batch_size==0)
         break;
       for(int i=0;i<curr_batch_size;i++)
         {   
            std::string read_seq(read_batch[i].read_seq);
            total_reads++;
            total_bases+=read_seq.length();
            std::vector<std::string> reads_without_ns;
            get_reads_spiliting_around_Ns(read_seq,read_seq.length(),reads_without_ns);
            for(int j=0;j<reads_without_ns.size();j++)
               {
                   
                  int nb_kmers_a=0;
                  pass_one_bloom_a(reads_without_ns[j],pass_one_arg->sample_kmers,nb_kmers_a);
                  total_kmers_a=total_kmers_a+nb_kmers_a;

                }
                             
                        
                                          
        } //for batch


  }//while
 pthread_mutex_lock( pass_one_arg->lock_S ) ;
 pass_one_arg->total_sample_kmers=pass_one_arg->total_sample_kmers+total_kmers_a;
 pass_one_arg->total_bases=pass_one_arg->total_bases+total_bases;
 pass_one_arg->total_reads=pass_one_arg->total_reads+total_reads;
 pthread_mutex_unlock( pass_one_arg->lock_S );
 free( read_batch ) ;
 pthread_exit(NULL);
 return NULL ;

}


void record_kmers_occured_bloom_a(std::vector<kmercode_length> kmers_list,bloom_util *bloomA,std::vector<bool> &occur,int &occured_kmers)
{

  int i;
  kmercode_length kmercode=0; 
  for(i=0;i<kmers_list.size();++i)
  {
     kmercode = kmers_list[i];
     if(bloomA->contain(kmercode))
     {             
        occur[i] = true ;
        occured_kmers++;

     }
     else
     {
        occur[i] = false ;
     }
  }

}

void record_trusted_read_positions(int read_length,std::vector<bool> occur,std::vector<int> threshold,std::vector<bool> &trusted_positions)
{
    int occurcnt = read_length- kmer_size + 1 ;
    int zerocnt = 0, onecnt = 0 ;

    for (int i = 0 ; i < read_length ; ++i )
    {
        if ( i >= kmer_size)
        {
            if ( occur[i - kmer_size] )
                --onecnt ;
            else
                --zerocnt ;
        }

        if ( i < occurcnt )
        {
            if ( occur[i] )
                ++onecnt ;
            else
                ++zerocnt ;
        }

        int sum = onecnt + zerocnt ;
        int adjust = 0 ;
        if ( onecnt > threshold[sum] + adjust )
        {

            trusted_positions[i] = true ;
        }
        else
        {
            trusted_positions[i] = false ;
        }

    }

}

void pass_two_bloom_b(std::vector<kmercode_length> kmers_list,std::vector<bool> trusted_positions,bloom_util *bloomB,
                      binary_store *solid_kmers,int &nb_total_kmers_b)
{
    
    int i=0, onecnt = 0;
    kmercode_length kmercode=0;
    for(i=0; i<kmer_size; ++i)
    {
        if ( trusted_positions[i] )
            ++onecnt ;
       
        if(onecnt == kmer_size)
        {
            kmercode= kmers_list[0];
            
           /* bloomB->add(kmercode);
            nb_total_kmers_b++;
            solid_kmers->write(&kmercode,sizeof(kmercode));*/
           
            if(bloomB->add(kmercode,true)==1)
            {
                 nb_total_kmers_b++;
                 solid_kmers->write(&kmercode,sizeof(kmercode));
            }
            
            

        }
    }


    for(i=1; i<kmers_list.size(); ++i)
    {
        if (trusted_positions[i+(kmer_size-1)])
            ++onecnt ;

        if ( trusted_positions[i - 1] )
            --onecnt ;

        if ( onecnt == kmer_size)
        {
            kmercode= kmers_list[i];
            
           /* bloomB->add(kmercode);
            nb_total_kmers_b++;
            solid_kmers->write(&kmercode,sizeof(kmercode));*/

            if(bloomB->add(kmercode,true)==1)
            {
                 nb_total_kmers_b++;
                 solid_kmers->write(&kmercode,sizeof(kmercode));
            }
            

        }
       
    }


}

void *pass_two_bloom_b_multi_threads( void *arg )
{

  struct threads_arg_store_trusted_kmers* pass_two_arg = (struct threads_arg_store_trusted_kmers*)arg;
  int max_batch_size= READ_BUFFER_PER_THREAD ;
  int curr_batch_size; 
  one_read *read_batch= (one_read*)malloc(sizeof(one_read) * max_batch_size) ;
  uint64_t total_kmers_b=0;uint64_t occured_kmers=0;
  while(1)
  {

       pthread_mutex_lock( pass_two_arg->lock_F ) ;
       curr_batch_size = pass_two_arg->reads->get_batch_of_reads(read_batch, max_batch_size);
       pthread_mutex_unlock( pass_two_arg->lock_F );
       if(curr_batch_size==0)
         break;
       for(int i=0;i<curr_batch_size;i++)
         {   
            std::string read_seq(read_batch[i].read_seq);
            std::vector<std::string> reads_without_ns;
            get_reads_spiliting_around_Ns(read_seq,read_seq.length(),reads_without_ns);
            for(int j=0;j<reads_without_ns.size();j++)
               {
                   int nb_kmers_b=0;
                   int nb_kmers_o=0;
                   int nb_kmers=reads_without_ns[j].length()-kmer_size+1;
                   std::vector<bool> occur(reads_without_ns[j].length()) ;
                   std::vector<bool> trusted_positions(reads_without_ns[j].length()) ;
                   std::vector<kmercode_length> kmers_list(nb_kmers);
                   get_kmers_one_read(reads_without_ns[j],kmers_list);
                   record_kmers_occured_bloom_a(kmers_list,pass_two_arg->sample_kmers,occur,nb_kmers_o);
                   record_trusted_read_positions(reads_without_ns[j].length(),occur,pass_two_arg->threshold,trusted_positions);
                   pass_two_bloom_b(kmers_list,trusted_positions,pass_two_arg->trusted_kmers,pass_two_arg->solid_kmers,nb_kmers_b);
                   total_kmers_b=total_kmers_b+nb_kmers_b;
                   occured_kmers=occured_kmers+nb_kmers_o;

               }
                             
                        
                                          
        } //for batch


  }//while
 pthread_mutex_lock( pass_two_arg->lock_S ) ;
 pass_two_arg->total_kmers_b=pass_two_arg->total_kmers_b+total_kmers_b;
 pass_two_arg->occured_kmers=pass_two_arg->occured_kmers+occured_kmers;
 pthread_mutex_unlock( pass_two_arg->lock_S );
 free( read_batch ) ;
 pthread_exit(NULL);
 return NULL ;


}//method


void *compute_branching_nodes_thread(void *arg)
{
  struct threads_arg_branching_kmers* branching_arg = (struct threads_arg_branching_kmers*)arg;
  int max_batch_size= READ_BUFFER_PER_THREAD ;
  int curr_batch_size; 
  kmercode_length *kmers_batch= (kmercode_length*)malloc(sizeof(kmercode_length) * max_batch_size) ;
  uint64_t total_kmers_branching=0;uint64_t total_kmers_file=0;
  while(1)
   {
       pthread_mutex_lock(branching_arg->lock_F);
       curr_batch_size = branching_arg->kmers_branching_obj->get_batch_of_kmers(kmers_batch, max_batch_size);
       pthread_mutex_unlock(branching_arg->lock_F);
       if(curr_batch_size==0)
          break;
       for(int i=0;i<curr_batch_size;i++)
          {
              kmercode_length kmercode=kmers_batch[i];
              total_kmers_file++;
             if(branching_arg->kmers_branching_obj->is_branching(kmercode))
                {
                    pthread_mutex_lock(branching_arg->lock_H);

                    if(branching_arg->kmers_branching_obj->ht.find(kmercode) == branching_arg->kmers_branching_obj->ht.end() )
                        {
                           branching_arg->kmers_branching_obj-> ht[kmercode]=1;
                           total_kmers_branching++;
                        }

                   else
                         branching_arg->kmers_branching_obj->ht[kmercode]=branching_arg->kmers_branching_obj->ht[kmercode]+1;

                   pthread_mutex_unlock(branching_arg->lock_H);

               }


         }//for

        if(curr_batch_size<max_batch_size)
           break;


   }//while

 pthread_mutex_lock(branching_arg->lock_S) ;
 branching_arg->total_kmers_file=branching_arg->total_kmers_file+total_kmers_file;
 branching_arg->total_kmers_branching=branching_arg->total_kmers_branching+total_kmers_branching;
 pthread_mutex_unlock(branching_arg->lock_S);
 free( kmers_batch ) ;
 pthread_exit(NULL);
 return NULL ;



}//method 

void help_me(int computed_coverage, double error_rate)
{

       std::cout<<"--- LightAssembler can not assemble your dataset !!! "<<std::endl;

       int start_gap_size=0;
       
       if (computed_coverage > 0)
        
       {

             if(computed_coverage >=280)
               { 
                      if(error_rate <=0.01)
                         start_gap_size=25;
                      else
                         start_gap_size=33;

               }

             else

             if(computed_coverage >=140)
               {

                      if(error_rate <=0.01)
                         start_gap_size=15;
                      else
                         start_gap_size=20;

               }

             else

             if(computed_coverage >=75)
              {

                      if(error_rate <=0.01)
                         start_gap_size=8;
                      else
                         start_gap_size=15;

              }
             else

              if(computed_coverage >=35)
              {

                      if(error_rate <=0.01)
                         start_gap_size=4;
                      else
                         start_gap_size=8;

              }

             else
              {
                     if(error_rate <=0.01)
                         start_gap_size=3;
                      else
                         start_gap_size=6;

              }


       }

       if((computed_coverage > 0)&&(gap_size < start_gap_size))
       std::cout<<"--- choose larger gap size, start with gap size g = "<<start_gap_size<<std::endl;

       if((computed_coverage <= 0))

       std::cout<<"--- errors in reading your dataset, average coverage = "<< computed_coverage <<std::endl;

     
       std::cout<<"--- maximum supported read length for this version = "<<MAX_READ_LENGTH<<std::endl;


       std::cout<<"--- try different values for k [kmer size] & g [gap size] or different dataset"<<std::endl;

       std::cout<<std::endl;
       exit(1);


}

void assemble(branching_kmers* kmers_branching_obj, uint64_t genomesize)
{
    uint64_t assembly_size=0;uint64_t number_contigs=0;uint64_t max_contig_length=0;
    kmers_branching_obj->clear_kmer_abundnce();

    kmers_branching_obj->start_iterator();
            
    traverse_kmers traversal(kmers_branching_obj->bloomB,kmers_branching_obj,10000000,500,20);
    //std::cout<<"------------------------------------------------------------------"<<std::endl;
    std::cout<<std::endl;
    //std::cout<<"--------------------------{ Graph traversal }---------------------"<<std::endl;
    std::cout<<"--- Graph traversal. "<<std::endl;
    std::cout<<std::endl;
    start_time(6);
    traversal.init(assembly_size, number_contigs, max_contig_length);
    end_time(6);
    double genome_coverage= static_cast<double>((static_cast<double>(assembly_size)/static_cast<double>(genomesize))*100.0);
    std::cout<<"--- number of contigs     = "<<number_contigs<<std::endl;
    std::cout<<"--- maximum contig length = "<<max_contig_length<<std::endl;
    std::cout<<"--- assembly size         = "<<assembly_size<<std::endl;
    std::cout<<"--- genome coverage       = "<<genome_coverage<< "%"<<std::endl;
    kmers_branching_obj->ht.clear();




}
void init(const std::vector<std::string> read_files,uint64_t genomesize,double error_rate, bool compute_g,int nb_threads,bool verbose)
{
   
    bloom_util bloomA((uint64_t)genomesize*1.5,0.01);
    bloom_util bloomB((uint64_t)genomesize*1.5,0.0005);
    pthread_attr_t pthread_attr ;
    pthread_t *threads = NULL;
    pthread_mutex_t mutex_sample_kmers ;
    pthread_mutex_t mutex_trusted_kmers_file ;
    pthread_mutex_t mutex_hashtable ;
    pthread_mutex_t mutex_sum ;
    
    if(nb_threads > 1)
      {
           pthread_attr_init( &pthread_attr ) ;
	   pthread_attr_setdetachstate( &pthread_attr, PTHREAD_CREATE_JOINABLE ) ;
	   threads = ( pthread_t * )malloc( sizeof( pthread_t ) * nb_threads ) ;
           pthread_mutex_init( &mutex_sample_kmers, NULL ) ;
           pthread_mutex_init( &mutex_trusted_kmers_file, NULL ) ;
           pthread_mutex_init( &mutex_sum, NULL ) ;
           pthread_mutex_init( &mutex_hashtable, NULL ) ;
           bloomA.set_nb_threads(nb_threads);
           bloomB.set_nb_threads(nb_threads);
     }

    uint64_t total_bases=0;uint64_t total_reads=0;
    std::string read_seq="";int read_length=0;
    double gapped_kmer=0.0;
    uint64_t  average_len=0;
    int computed_coverage=0;

    reads_parsing reads(read_files);
   
   if(compute_g||(gap_size<=0))
    {
       //std::cout<<"------------------------------------------------------------------"<<std::endl;
       //std::cout<<"-------------------{ Parameters extrapolation }-------------------"<<std::endl;
       std::cout<<"--- Parameters extrapolation. "<<std::endl;
       std::cout<<std::endl;
       start_time(1);
       if(nb_threads==1)
       {
         
           while(reads.get_next_sequence(read_seq,read_length))
           {

                  total_reads++;
                  total_bases+=read_seq.length();

           }

       }
      else
       {

            struct threads_arg_parameters_extrapolation arg;
            void *pthread_status ;
            arg.total_bases=0;
            arg.total_reads=0;
            arg.reads=&reads;
            arg.lock_A=&mutex_sample_kmers;
            arg.lock_S=&mutex_sum;

            for ( int i = 0 ; i < nb_threads ; ++i )
            {
              pthread_create( &threads[i], &pthread_attr, pass_zero_parameters_extrapolation_multi_threads, (void *)&arg ) ;

            }
            for ( int i = 0 ; i < nb_threads ; ++i )
            {
               pthread_join( threads[i], &pthread_status ) ;
            }

           total_bases=arg.total_bases;
           total_reads=arg.total_reads;

       }
       end_time(1);

       if((total_reads >0))
       {

          average_len=total_bases/total_reads;
          if((total_bases%total_reads)>(total_reads/2))
              average_len++;

          computed_coverage=(total_bases/genomesize);


          if(computed_coverage >=280)
             {
                  if(error_rate <=0.01)
                         gap_size=25;
                  else
                         gap_size=33;

             }

          else
 
          if(computed_coverage >=140)
             {

                   if(error_rate <=0.01)
                          gap_size=15;
                   else
                          gap_size=20;

             }

          else

          if(computed_coverage >=75)
             {

                   if(error_rate <=0.01)
                          gap_size=8;
                   else
                          gap_size=15;

             }
          else

          if(computed_coverage >=35)
            {

                   if(error_rate <=0.01)
                          gap_size=4;
                   else
                          gap_size=8;

            }

          else
            {
                   if(error_rate <=0.01)
                          gap_size=3;
                   else
                          gap_size=6;

            }

           gapped_kmer= static_cast<double>(1.0/static_cast<double>(gap_size));


        }
       else
        
          help_me(computed_coverage,error_rate);


       if(verbose)
       {                
                  std::cout<<"--- start with gap size g = "<<gap_size<< std::endl;

                  //std::cout<<"--- probability of a kmer added to BloomA = "<<gapped_kmer<< std::endl;

                  std::cout<<"--- average read length = "<<average_len<< std::endl;

                  std::cout<<"--- average sequencing coverage = "<<computed_coverage<<std::endl;

       }

                 //  std::cout<<"------------------------------------------------------------------"<<std::endl;

                  std::cout<<std::endl;


   }

    reads.rewind_all();
    uint64_t total_kmers_samples=0;
    //std::cout<<"------------------------------------------------------------------"<<std::endl;
    //std::cout<<"-------------------{ Uniform kmers sampling }---------------------"<<std::endl;
    std::cout<<"--- Uniform kmers sampling. "<<std::endl;
    std::cout<<std::endl;
    start_time(2);
    if(nb_threads==1)
    {
           while(reads.get_next_sequence(read_seq,read_length))
           {        

                  total_reads++;
                  total_bases+=read_seq.length();
                  std::vector<std::string> reads_without_ns;
                  get_reads_spiliting_around_Ns(read_seq,read_length,reads_without_ns);
                  for(int i=0;i<reads_without_ns.size();i++)
                  {
                      int nb_samples_r=0;
                      pass_one_bloom_a(reads_without_ns[i],&bloomA,nb_samples_r);
                      total_kmers_samples=total_kmers_samples+nb_samples_r;
                  }


            }
    }

    else
    {
      struct threads_arg_sample_kmers arg;
      void *pthread_status ;
      arg.total_sample_kmers=0;
      arg.total_bases=0;
      arg.total_reads=0;
      arg.sample_kmers=&bloomA;
      arg.reads=&reads;
      arg.lock_A=&mutex_sample_kmers;
      arg.lock_S=&mutex_sum;

     for ( int i = 0 ; i < nb_threads ; ++i )
	{
	      pthread_create( &threads[i], &pthread_attr, pass_one_bloom_a_multi_threads, (void *)&arg ) ;
            
	}
     for ( int i = 0 ; i < nb_threads ; ++i )
	{
		
              pthread_join( threads[i], &pthread_status ) ;
              
	}

     total_kmers_samples=arg.total_sample_kmers;
     total_bases=arg.total_bases;
     total_reads=arg.total_reads;


    }
    end_time(2);
    double bloom_fpr= bloomA.get_false_positive_rate();
    std::vector<double> untrust(kmer_size+1);
    std::vector<int> threshold(kmer_size+1);
    if((!compute_g)&&(total_reads >0))
    {
              gapped_kmer= static_cast<double>(1.0/static_cast<double>(gap_size));
              average_len=total_bases/total_reads;
              if((total_bases%total_reads)>(total_reads/2))
              average_len++;
              computed_coverage=(total_bases/genomesize);
    }

   if(computed_coverage <=0)

       help_me(computed_coverage,error_rate);

   if(verbose)
    {
         std::cout<<"--- total number of kmers in BloomA = "<<total_kmers_samples <<std::endl;
    
         std::cout<<"--- BloomA false positive rate = "<<bloom_fpr<< std::endl; 
        
         if(!compute_g)
           {
    
                 //std::cout<<"--- probability of a kmer added to BloomA = "<<gapped_kmer<< std::endl;
    
                 std::cout<<"--- average read length = "<<average_len<< std::endl;

                 std::cout<<"--- average sequencing coverage = "<<computed_coverage<<std::endl;

           }

    }
   
   if(bloom_fpr > 0.65)
   {
        
       help_me(computed_coverage,error_rate);

   }

   compute_distribution_untrusted_k_positions(untrust,computed_coverage,error_rate,gapped_kmer,bloom_fpr,verbose);
   compute_threshold_k_positions(untrust,threshold,gapped_kmer);
   uint64_t total_kmers_b=0;uint64_t occured_kmers=0;
   binary_store solid_kmers(return_file_name(solid_kmers_file),sizeof(kmercode_length),true);
   std::cout<<std::endl;
   //std::cout<<"------------------------------------------------------------------"<<std::endl;
   //std::cout<<"---------------{ Trusted/untrusted kmers filtering }--------------"<<std::endl;
   std::cout<<"--- Trusted/untrusted kmers filtering. "<<std::endl;
   std::cout<<std::endl;
   reads.rewind_all();
   start_time(4);
   if(nb_threads==1)
   {
      read_seq="";read_length=0;
      while(reads.get_next_sequence(read_seq,read_length))
        {
            std::vector<std::string> reads_without_ns;
            get_reads_spiliting_around_Ns(read_seq,read_length,reads_without_ns);
            for(int i=0;i<reads_without_ns.size();i++)
            {
                   
               int nb_kmers_b=0;
               int nb_kmers_o=0;
               int nb_kmers=reads_without_ns[i].length()-kmer_size+1;
               std::vector<bool> occur(reads_without_ns[i].length()) ;
               std::vector<bool> trusted_positions(reads_without_ns[i].length()) ;
               std::vector<kmercode_length> kmers_list(nb_kmers);
               get_kmers_one_read(reads_without_ns[i],kmers_list);
               record_kmers_occured_bloom_a(kmers_list,&bloomA,occur,nb_kmers_o);
               record_trusted_read_positions(reads_without_ns[i].length(),occur,threshold,trusted_positions);
               pass_two_bloom_b(kmers_list,trusted_positions,&bloomB,&solid_kmers,nb_kmers_b);
               total_kmers_b=total_kmers_b+nb_kmers_b;
               occured_kmers=occured_kmers+nb_kmers_o;
            }

        }

   solid_kmers.close();
   reads.close();
  }
   else
   {
    struct threads_arg_store_trusted_kmers arg;
    void *pthread_status ;
    arg.threshold=threshold;
    arg.total_kmers_b=0;
    arg.occured_kmers=0;
    arg.sample_kmers=&bloomA;
    arg.trusted_kmers=&bloomB;
    arg.reads=&reads;
    arg.solid_kmers=&solid_kmers;
    arg.lock_F=&mutex_trusted_kmers_file;
    arg.lock_S=&mutex_sum;
    
    for ( int i = 0 ; i < nb_threads ; ++i )
	{
	      pthread_create( &threads[i], &pthread_attr, pass_two_bloom_b_multi_threads, (void *)&arg ) ;
            
	}
    for ( int i = 0 ; i < nb_threads ; ++i )
	{
		
              pthread_join( threads[i], &pthread_status ) ;
              
	}
   
    
   total_kmers_b=arg.total_kmers_b;
   occured_kmers=arg.occured_kmers;

   solid_kmers.close();
   reads.close();

  }  

  end_time(4);

  if(verbose)
  {
    //std::cout<<"--- total number of kmers occured in BloomA + false positive "<< occured_kmers <<std::endl;
    std::cout<<"--- total number of kmers in BloomB = "<<total_kmers_b <<std::endl;
    std::cout<<"--- BloomB false positive rate = "<<bloomB.get_false_positive_rate()<< std::endl;

  }

  if((total_kmers_b > 0)&&(bloomB.get_false_positive_rate()< .009))// 0.001) //.09)
     {
       //std::cout<<"------------------------------------------------------------------"<<std::endl;
       std::cout<<std::endl;
       //std::cout<<"------------------{ Branching-kmers computation }-----------------"<<std::endl;
       std::cout<<"--- Branching-kmers computation. "<<std::endl;
       std::cout<<std::endl;
       uint64_t total_kmers_branching=0;
       kmercode_length kmercode=0;
       binary_store solid_kmers_t(return_file_name(solid_kmers_file),sizeof(kmercode),false);
       hash_table ht;
       branching_kmers* kmers_branching_obj=new branching_kmers(solid_kmers_t,bloomB,ht);
       start_time(5);

        if(nb_threads==1)
           kmers_branching_obj->compute_branching_kmers(total_kmers_branching);
        else
          {

            struct threads_arg_branching_kmers arg;
            void *pthread_status ;
            arg.total_kmers_branching=0;
            arg.total_kmers_file=0;
            arg.kmers_branching_obj=kmers_branching_obj;
            arg.lock_H=&mutex_hashtable;
            arg.lock_F=&mutex_trusted_kmers_file;
            arg.lock_S=&mutex_sum;

            for ( int i = 0 ; i < nb_threads ; ++i )
	       {
                  pthread_create( &threads[i], &pthread_attr, compute_branching_nodes_thread, (void *)&arg ) ;
            
	       }
            for ( int i = 0 ; i < nb_threads ; ++i )
	       {
		
                     pthread_join( threads[i], &pthread_status) ;
              
	       }

           

           /* if(verbose)
                   std::cout<<"--- total number of solid kmers in the file: "<< arg.total_kmers_file <<std::endl;*/

             total_kmers_branching=arg.total_kmers_branching;

           }

           end_time(5);
           if(total_kmers_branching==0)
              {

                   std::cerr <<"--- no branching kmers found in this data set."<<std::endl;
                   solid_kmers_t.close();
                   delete_created_file(return_file_name(solid_kmers_file));
                   help_me(computed_coverage,error_rate);

              }
           else
              {
                   if(verbose)
                   std::cout<<"--- number of branching kmers = "<<total_kmers_branching<<std::endl;
              }

           solid_kmers_t.close();
           delete_created_file(return_file_name(solid_kmers_file));
           assemble(kmers_branching_obj, genomesize);


     }
  else
    {


       delete_created_file(return_file_name(solid_kmers_file));
       help_me(computed_coverage,error_rate);


    }     

}
#endif
