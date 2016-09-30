#include "ReadsParsing.hpp"
reads_parsing::reads_parsing(std::vector<std::string> file_names)
{
     nb_files=file_names.size();
     if( nb_files >= MAX_READ_FILE )
         {
            std::cerr<<"--- the number of read files exceeds the limit "<<MAX_READ_FILE<<std::endl;
            exit( 1 ) ;
         }

        read_files = (read_file**) malloc(sizeof(read_file *)*nb_files);
        
         for(int i=0;i<nb_files;i++)
         {

            const char * f_name = file_names[i].c_str();
            gzFile fp = gzopen(f_name,"r");
            if(fp == 0)
            {
              std::cerr<<"--- can not open the file "<<file_names[i]<<std::endl;
              exit(1);
            }
            char delim =(char)gzgetc(fp);
            if(delim=='>' || delim=='@')
            { // file in a fasta/q format
              gzclose(fp);
            }
            else
            {
              std::cerr<<"--- can not recognize the file format "<<file_names[i]<<std::endl;
              exit(1);

            }

            read_files[i] = (read_file *)calloc(1, sizeof(read_file));
            read_files[i]->file_name=f_name;
         }
         rewind_all();
         file_names.clear();


}
reads_parsing::~reads_parsing()
{
   for(int i=0;i<nb_files;i++)
       {
         free(read_files[i]);
       }

}

void reads_parsing::open_file(int file_indx)
{
     read_files[file_indx]->fp = gzopen(read_files[file_indx]->file_name,"r");
     read_files[file_indx]->seq= kseq_init(read_files[file_indx]->fp);

}
void reads_parsing::close_file(int file_indx)
{
     kseq_destroy(read_files[file_indx]->seq);
     gzclose(read_files[file_indx]->fp);
     read_files[file_indx]->fp= NULL;


 }
void reads_parsing::close()
 {
     for(int i=0;i<nb_files;i++)
       {
          if(read_files[i]->fp !=NULL)
           {
             kseq_destroy(read_files[i]->seq);
             gzclose(read_files[i]->fp);
             read_files[i]->fp =NULL;
           }
           //free(read_files[i]);
       }
   }



void reads_parsing::rewind_all()
 {
      for(int i=0;i<nb_files;i++)
       {
          if(read_files[i]->fp !=NULL)
           {
             kseq_destroy(read_files[i]->seq);
             gzclose(read_files[i]->fp);
             read_files[i]->fp =NULL;
           }
       }

     current_file_indx = 0;
     open_file(current_file_indx);
 }


 bool reads_parsing::get_next_seq_from_file(std::string &seq_read, int &len,int file_indx)
 {
      int l=kseq_read(read_files[file_indx]->seq);
      if(l<0)
       return false;

      if(read_files[file_indx]->seq->seq.l > MAX_READ_LENGTH)
       {
          //std::cerr<<"--- maximum supported read length for this version = "<<MAX_READ_LENGTH<<std::endl;
            return false;
       }

      seq_read=read_files[file_indx]->seq->seq.s;
           len=read_files[file_indx]->seq->seq.l;
      return true;

  }
 bool reads_parsing::get_next_seq_from_file_buffer(char *seq_read, int &len,int file_indx)
 {
      int l=kseq_read(read_files[file_indx]->seq);
      if(l<0)
       return false;

      if(read_files[file_indx]->seq->seq.l > MAX_READ_LENGTH)
       {       
         // std::cerr<<"--- maximum supported read length for this version = "<<MAX_READ_LENGTH<<std::endl;
            return false;
       }

      strcpy(seq_read,read_files[file_indx]->seq->seq.s);
      len=read_files[file_indx]->seq->seq.l;
      return true;

  }
bool reads_parsing::get_next_sequence(std::string &seq_read, int &len)
 {
       bool success = get_next_seq_from_file(seq_read,len,current_file_indx);
       if (success)
        return true;
    
       if ( current_file_indx < nb_files-1 )
       {
            close_file(current_file_indx);
            current_file_indx++;
            open_file(current_file_indx);
            return get_next_sequence(seq_read,len);
       }
       return false;
 }

int reads_parsing::get_batch_of_reads(one_read *reads,int max_batch_size)
{
  int batch_size=0;
  while (batch_size < max_batch_size)
   {
      int read_len=0;
      bool success = get_next_seq_from_file_buffer(reads[batch_size].read_seq,read_len,current_file_indx);
      if (success)
        { 
           batch_size++;
           continue;

        }


      if ( current_file_indx < nb_files-1 )
       {
            close_file(current_file_indx);
            current_file_indx++;
            open_file(current_file_indx);
            continue;
       }
       
       if((!success)&& batch_size>0)
       return batch_size;
       else
       return 0;
   }
   return batch_size;
}
