#ifndef READSPARSING_H
#define READSPARSING_H
#include <iostream>
#include <zlib.h>
#include "kseq.h"
#include "KmerUtil.hpp"
#include "Utility.hpp"
#include <vector>
#include <string>

static const int MAX_READ_FILE = 100;
KSEQ_INIT(gzFile, gzread)
struct read_file
{
 gzFile fp;
 kseq_t *seq;
 const char *file_name;
};
struct one_read
{
 char read_seq[MAX_READ_LENGTH] ;
};

class reads_parsing
{
   
   public:
   read_file **read_files;
   int nb_files; 
   int current_file_indx ;
   reads_parsing(std::vector<std::string> file_names);
  ~reads_parsing();
   void open_file(int file_indx);
   void close_file(int file_indx);
   void close();
   void rewind_all();
   bool get_next_seq_from_file(std::string &seq_read, int &len,int file_indx);
   bool get_next_seq_from_file_buffer(char *seq_read, int &len,int file_indx);
   bool get_next_sequence(std::string &seq_read, int &len);
   int  get_batch_of_reads(one_read *reads,int max_batch_size);
};

#endif
