#include <iostream>
#include "KmerUtil.hpp"
#include "KmersScanning.hpp"
#include <vector>
#include <getopt.h>//getopt_long
#include <unistd.h>
#include <sys/stat.h>//stat
#include <pthread.h>
#include <sstream>
struct program_options
{
    int  k;
    int  g;
    double e;
    int  t;
    uint64_t genomesize;
    std::vector<std::string> read_files;
    std::string output_prefix_name;
    bool verbose;
    bool compute_g;
    program_options():k(31),g(0),e(0.01),t(1),compute_g(true),
                                              genomesize(0),output_prefix_name("LightAssembler"),verbose(false){}
};
void print_usage()
{
    std::cerr <<std::endl;
    std::cerr <<" ********************************** <<<   LightAssembler  >>> *************************************** "<<std::endl<<std::endl ;
    std::cerr <<" Light Version of an assembly algorithm for short reads in FASTA/FASTQ/FASTA.gz/FASTQ.gz formats."<<std::endl<<std::endl ;
    std::cerr <<" Usage: ./LightAssembler [Options] ...FASTA/FASTQ/FASTA.gz/FASTQ.gz files"<<std::endl;
    std::cerr <<std::endl<<
    "  [-k] kmer size                                [default: 31] "<<std::endl<<
    "  [-g] gap size                                 [default: 25X:3 35X:4 75X:8 140X:15 280X:25] "<<std::endl<<
    "  [-e] expected error rate                      [default: 0.01 ] "<<std::endl<<
    "  [-G] genome size                              [default: 0] "<<std::endl<<
    "  [-t] number of threads                        [default: 1]"<<std::endl<<
    "  [-o] output file name                         [default: LightAssembler] "<<std::endl;
    std::cerr <<std::endl;
    std::cerr <<" Typical LightAssembler Command Line :"<<std::endl<<std::endl;
    std::cerr <<" ./LightAssembler -k 31 -g 13 -e 0.01 -G 1000000 -t 1 -o LightAssembler read_file1 read_file2 --verbose "<<std::endl;
    std::cerr <<std::endl;
    std::cerr <<" **************************************************************************************************** "<<std::endl<<std::endl ;

}
void parse_options(int argc,char **argv,program_options &opt)
{
    int verbose_flag=0;
    const char* opt_string =":k:g:e:G:o:t:";
    static struct option long_options[]=
    {
        {"verbose",no_argument,&verbose_flag,1},
        {"kmer size",required_argument,NULL,'k'},
        {"gap size",required_argument,NULL,'g'},
        {"error rate",required_argument,NULL,'e'},
        {"threads",required_argument,NULL,'t'},
        {"genome size",required_argument,NULL,'G'},
        {"output file name",required_argument,NULL,'o'},
        {NULL,no_argument,NULL,0}
    };
   int option_index=0;
   int ch;
   while((ch=getopt_long(argc,argv,opt_string,long_options,&option_index))!=-1)
   {
        if(ch=='k')
        {
            std::istringstream argument(optarg);
            if(!(argument>>opt.k)||!(argument.eof()))
             {
                std::cerr<<std::endl;
                std::cerr<<"--- invalid argument of -k "<<optarg<<std::endl;
                print_usage();
                exit(1);
             }

        }       
        else if(ch=='g')
        {
            std::istringstream argument(optarg);
            if(!(argument>>opt.g)||!(argument.eof()))
             {
                std::cerr<<std::endl;
                std::cerr<<"--- invalid argument of -g "<<optarg<<std::endl;
                print_usage();
                exit(1);
             }
            opt.compute_g=false;
        }
        else if(ch=='e')
        {
            std::istringstream argument(optarg);
            if(!(argument>>opt.e)||!(argument.eof()))
             {
                std::cerr<<std::endl;
                std::cerr<<"--- invalid argument of -e "<<optarg<<std::endl;
                print_usage();
                exit(1);
             }

        }
        else if(ch=='t')
        {
            std::istringstream argument(optarg);
            if(!(argument>>opt.t)||!(argument.eof()))
             {
                std::cerr<<std::endl;
                std::cerr<<"--- invalid argument of -t "<<optarg<<std::endl;
                print_usage();
                exit(1);
             }

        }
        else if(ch=='G')
        {

            std::istringstream argument(optarg);
            if(!(argument>>opt.genomesize)||!(argument.eof()))
             {
                std::cerr<<std::endl;
                std::cerr<<"--- invalid argument of -G "<<optarg<<std::endl;
                print_usage();
                exit(1);
             }

       }
       else if(ch=='o')
       {
           opt.output_prefix_name=optarg;
          
       }
       else if(ch=='?')
       {

             std::cerr<<std::endl;
             if (isprint (optopt))
             std::cerr<<"--- invalid option -"<< static_cast<char>(optopt) <<std::endl;
             else
             std::cerr<<"--- invalid option "<< argv[optind-1]<<std::endl;
             print_usage();
             exit(1);
       }
       else if(ch==':')
       {

             std::cerr<<std::endl;
             std::cerr<<"--- missing argument of -"<< static_cast<char>(optopt)<<std::endl;
             print_usage();
             exit(1);

       }

  }

  for(int i=optind;i<argc;++i)
  {opt.read_files.push_back(argv[i]);}
  if(verbose_flag)
  { opt.verbose=true;}

}
bool check_options(program_options &opt)
{
    bool success=true;
    std::cerr<<std::endl;
    //************************** kmer size**************************************************
    if(opt.k>0)
    {
       if(opt.k>((int)sizeof(kmercode_length)*4))
          {
          
           std::cerr<<"--- maximum support kmer size for this compiled version : "<<(sizeof(kmercode_length)*4)<<std::endl;
           std::cerr<<"--- use make k="<<opt.k<<std::endl;
           success=false;
         
          }
       else if(opt.k%2==0)
          {opt.k--;std::cout<<"--- to avoid palindromes, kmer size must be odd, suggested kmer size "<<opt.k<<std::endl;}
    }
    else
    {std::cerr<<"--- invalid value for kmer size: "<<opt.k<<std::endl;success=false;}
    //************************** gap size**************************************************
   /* if(opt.g>0)
    {
        if(opt.g>opt.k)
            {std::cerr<<"--- supported gap sizes 1 ~ kmer_size: "<<std::endl;success=false;}
    }
    else */
   /*
    if(opt.g<=0)
    {std::cerr<<"--- invalid value for gap size: "<<opt.g<<std::endl;success=false;}*/
    //************************** coverage **************************************************
    if(opt.e>=1)
    {std::cerr<<"--- invalid value for expected error rate : "<<opt.e<<std::endl;success=false;}
    //************************** threads **************************************************
    if(opt.t<=0)
    {std::cerr<<"--- invalid value for number of threads "<<opt.t<<" , value must be at least 1"<<std::endl;success=false;}

    //************************** read files**************************************************
    if(opt.read_files.size()==0)
    {std::cerr<<"--- no read files specified as inputs"<<std::endl;success=false;}
    else
    {
        struct stat stat_file_info;
        int int_stat;
        std::vector<std::string>::const_iterator it;
        for(it=opt.read_files.begin();it != opt.read_files.end();++it)
        {
            int_stat=stat(it->c_str(),&stat_file_info);
            if(int_stat != 0)
            {
             std::cerr<<"--- error: file not found "<<*it<<std::endl;
             success=false;

            }//if
        }//for


    }//else

//************************** genome size **********************************************************
    if(opt.genomesize<=0)
    {std::cerr<<"--- invalid value for genome size: "<<opt.genomesize<<std::endl;success=false;}

   return success;

}

int main(int argc,char** argv)
{
    program_options opt;
    parse_options(argc,argv,opt);
    if(!check_options(opt))
    {
     print_usage();
     exit(1);
    }
    kmer_size=opt.k;
    gap_size=opt.g;
    kmer_mask=((static_cast<kmercode_length>(1))<<(kmer_size*2))-1;
    prefix=opt.output_prefix_name;
    start_time(0);
    init(opt.read_files,opt.genomesize,opt.e,opt.compute_g,opt.t,opt.verbose);
    std::cout<<std::endl;
    //std::cout<<"------------------------------------------------------------------"<<std::endl;
    //std::cout<<"-------------------{ Assembly process finished }------------------"<<std::endl;
    std::cout<<"--- The assembly session is finished. "<<std::endl;
    std::cout<<std::endl;
    //std::cout<<"------------------------------------------------------------------"<<std::endl;
    end_time(0);
    //std::cout<<"------------------------------------------------------------------"<<std::endl;
    std::cout<<std::endl;
    return 0;

}
