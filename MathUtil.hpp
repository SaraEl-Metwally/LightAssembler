#ifndef MATHUTIL_H
#define MATHUTIL_H

void get_cumulative_binomial_distribution( std::vector<double> & F, int l, double p )
{
    // p is the probability of getting 1.
    int i ;
    double coef = 1 ;
    double exp = pow( 1 - p, l ) ;
    F[0] = pow( 1 - p, l ) ;
    
    for ( i = 1 ; i <= l ; ++i )
    {
        coef = coef / i * ( l - i + 1 ) ;
        exp =  exp / ( 1 - p ) * p ;
        F[i] = F[i - 1] + coef * exp ;
    }
    
}




void compute_distribution_untrusted_k_positions(std::vector< double> & untrust,int coverage,double error_rate, 
                                                                                double alpha, double false_positive, bool verbose)
{

    double prob_incorr_kmer_sample=(1-exp(-1*alpha*error_rate*coverage));
    double prob_incorr_kmer_sample_with_fb=(false_positive+((1-false_positive)*prob_incorr_kmer_sample));
    if(verbose)
    std::cout<<"--- probability of an incorrect kmer appears in the sample : "<<prob_incorr_kmer_sample_with_fb<<std::endl;
    get_cumulative_binomial_distribution(untrust, kmer_size,prob_incorr_kmer_sample_with_fb) ;

}

void compute_threshold_k_positions(std::vector<double> untrust, std::vector<int> &threshold, double alpha)
{
    int x = floor(kmer_size*alpha);
    int tk=0;
    for(int j=0;j<=kmer_size;j++)
    {

            if(untrust[j]>= 0.9)
             {  
                  
                    int tx= j+x;
                    if(tx >= kmer_size)
                        tx=j;

                     tk=tx;
                     threshold[kmer_size]=tk;

                     break;
             }

    }

    for(int i=1;i<=kmer_size;i++)
    {
             
           double t1  = static_cast<double>(static_cast<double>(i)/static_cast<double>(kmer_size));
           double t2  = (t1*tk);
           threshold[i]=round(t2);

    }

}

#endif
