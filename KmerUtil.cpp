#include "KmerUtil.hpp"
int kmer_size;
int gap_size;
kmercode_length kmer_mask;

int nt2int(char nt)
{
    
    if(nt=='A'||nt=='a')
        return 0;
    if(nt=='C'||nt=='c')
        return 1;
    if(nt=='G'||nt=='g')
        return 2;
    if(nt=='T'||nt=='t')
        return 3;
}
int rev_nt2int(int nt)
{
    int rev_nt= ~nt;
    return rev_nt &3;
}


int first_nt(kmercode_length kmer)
{

    int result;
    #ifdef largeintlib
    LargeInt<kmer_precision> code = kmer;
    result = code.toInt()&3;
    #else
    result = kmer&3;
    #endif
    return result;
}

int code2nt(kmercode_length code, int which_nt)
{
    kmercode_length temp = code;
    temp = temp >> (2*(kmer_size-1-which_nt));
    return first_nt(temp&3);
}

kmercode_length get_canonical_kmer_code(kmercode_length code)
{
    int i ;
    kmercode_length rev_com = get_reverse_complement(code);
    return rev_com < code ? rev_com : code ;
}
kmercode_length get_reverse_complement(kmercode_length code)
{
    int i ;
    kmercode_length rev_comp =  (static_cast<kmercode_length>(0)) ;
    for ( i = 0 ; i < kmer_size ; ++i )
    {
        kmercode_length tmp = ( code >> ( 2 * i ) ) & (static_cast<kmercode_length>(3)) ;
        rev_comp = ( rev_comp << 2 ) | ( (static_cast<kmercode_length>(3)) - tmp ) ;
    }


    return rev_comp;
}
int code2seq (kmercode_length code, char seq[])
{
    int i;
    kmercode_length temp = code;
    char bin2nt[4] = {'A','C','G','T'};

    for (i=kmer_size-1; i>=0; i--)
    {
        seq[i]=bin2nt[first_nt(temp&3)];
        temp = temp>>2;
    }
    seq[kmer_size]='\0';
    return kmer_size;
}
kmercode_length next_kmer(kmercode_length kmer_code, int added_nt, int *strand)
{
  

    kmercode_length new_kmer_code=0;
    kmercode_length temp_kmer_code=0;
    kmercode_length revcomp_new_kmer=0;

    if (strand != NULL && *strand == 1)
        temp_kmer_code = get_reverse_complement(kmer_code);
    else
        temp_kmer_code = kmer_code;

    new_kmer_code = (((temp_kmer_code) * 4 )  + added_nt) & kmer_mask;
    revcomp_new_kmer = get_reverse_complement(new_kmer_code);
    if (strand != NULL)
        *strand = (new_kmer_code < revcomp_new_kmer)?0:1;

    return std::min(new_kmer_code,revcomp_new_kmer);
}
void char_revcomp(char nt,char &cnt)
{
    if(nt=='A'||nt=='a')
        {cnt='T'; return;}
    if(nt=='C'||nt=='c')
        {cnt='G';return;}
    if(nt=='G'||nt=='g')
        {cnt='C';return;}
    if(nt=='T'||nt=='t')
        {cnt='A';return;}
}
void revcomp_sequence(std::vector<char> &sequence, int len)
{
		  int i;
		  unsigned char t;
		  for (i=0;i<len/2;i++)
		  {
			  t=sequence[i];
			  char_revcomp(sequence[len-i-1],sequence[i]);
			  char_revcomp(t,sequence[len-i-1]);
		  }
		  if (len%2==1)
			  char_revcomp(sequence[len/2],sequence[len/2]);

}
void get_kmers_one_read(std::string read, std::vector<kmercode_length> &kmers_list)
{

    int i;
    kmercode_length kmercode=0,kmercode_can=0;
    for(i=0; i<kmer_size; ++i)
    {

        kmercode=kmercode*4+nt2int(read[i]);


    }
    kmercode_can=get_canonical_kmer_code(kmercode);
    kmers_list[i - kmer_size]=kmercode_can;
    for(i=1; i<read.length()-kmer_size+1; ++i)
    {
        kmercode=(kmercode*4+nt2int(read[i+(kmer_size-1)])& kmer_mask);
        kmercode_can=get_canonical_kmer_code(kmercode);
        kmers_list[i]=kmercode_can;
    }

}

float needleman_wunch(std::string a, std::string b)
{
    float gap_score = -5;
    float mismatch_score = -5;
    float match_score = 10;
    #define nw_score(x,y) ( (x == y) ? match_score : mismatch_score )
    int n_a = a.length(), n_b = b.length();
    std::vector< std::vector<float> > score(n_a+1,std::vector<float>(n_b+1));

    for (int i = 0; i <= n_a; i++)
        score[i][0] = gap_score * i;
    for (int j = 0; j <= n_b; j++)
        score[0][j] = gap_score * j;


    for (int i = 1; i <= n_a; i++)
    {
        for (int j = 1; j <= n_b; j++)
        {
            float match = score[i - 1][j - 1] + nw_score(a[i-1],b[j-1]);
            float del =  score[i - 1][j] + gap_score;
            float insertion = score[i][j - 1] + gap_score;
            score[i][j] = std::max(std::max(match, del), insertion);
        }
    }

    // traceback
    int i=n_a, j=n_b;
    float similarity = 0;
    while (i > 0 && j > 0)
    {
        float score_current = score[i][j], score_diagonal = score[i-1][j-1], score_up = score[i][j-1], score_left = score[i-1][j];
        if (score_current == score_diagonal + nw_score(a[i-1], b[j-1]))
        {
            if (a[i-1]== b[j-1])
                similarity++;
            i -= 1;
            j -= 1;
        }
        else
        {
            if (score_current == score_left + gap_score)
                i -= 1;
            else if (score_current == score_up + gap_score)
                    j -= 1;
        }
    }
    similarity /= std::max( n_a, n_b);

    return similarity;
}

void get_reads_spiliting_around_Ns(std::string read_seq,int read_length,std::vector<std::string> & reads)
{
          int i=0;
          int indx=0;
 
            while (i < read_length)
            {
                indx=0;

                while (read_seq[i] =='N' && i < read_length)
                {
                    i++;
                }
                while ( (read_seq[i+indx] !='N') &&  ((i +indx) < read_length))
                {
                    indx++;
                }
                std::string read_without_ns = read_seq.substr(i,indx);
                if(read_without_ns.length()>=kmer_size)
                  {
                    reads.push_back(read_without_ns);
                  }
                i = indx+i;
            }
 
             

}
