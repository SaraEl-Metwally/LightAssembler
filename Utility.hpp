#ifndef UTILITY_H
#define UTILITY_H
#include <string>
#include <ctime>
#include <cstdio>
#include <iostream>
const   std::string log_file="assembly_log_details";
const   std::string solid_kmers_file="solid_kmers";
const   std::string assembly_file="contigs.fasta";
extern  std::string prefix;// In case if you asked the user about the prefix name for his files.
extern  time_t rawtime ;
static const int SECS_PER_MIN = 60 ;
static const int SECS_PER_HOUR = 3600;
static const int MAX_READ_LENGTH = 1024;
static const int READ_BUFFER_PER_THREAD = 1024;

#define start_time(t) \
time(&rawtime);\
double start_time ## t = rawtime;

#define end_time(t)\
time(&rawtime);\
double end_time  ## t  = rawtime;\
double elspsed   ## t  = difftime(end_time ## t,start_time ## t);\
int hours        ## t  = elspsed ## t / SECS_PER_HOUR;\
int minutes      ## t  = elspsed ## t / SECS_PER_MIN;\
int mins_left    ## t  = minutes ## t % SECS_PER_MIN;\
int secs_left    ## t  = (int)elspsed ## t % SECS_PER_MIN; \
std::cout<<"--- h("<<hours ## t<<")"<<":m("<<mins_left ## t<<")"<<":s("<<secs_left ## t<<") elapsed time."<< std::endl;

std::string return_file_name(const std::string file_name);
void delete_created_file(const std::string file_name);
#endif
