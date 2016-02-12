## LightAssembler
Lightweight resources assembly algorithm for high-throughput sequencing reads.
#### System requirements 
64-bit machine with g++ compiler or gcc in general.
#### Installation 
1. Clone the [GitHub repo](https://github.com/SaraEl-Metwally/LightAssembler), e.g. with `git clone https://github.com/SaraEl-Metwally/LightAssembler.git`
2. Run `make` in the repo directory for **k <= 31**  or `make k=kmersize` for **k > 31**.

#### Quick usage guide
``` ./LightAssembler -k [kmer size] -g [gap size] -e [error rate] -G [genome size] -t
[threads] -o [output prefix] [input files] --verbose ``` 

``` 
* [-k] kmer size                [default: 31]
* [-g] gap size                 [default: 25X:3 35X:4 75X:8 140X:15 280X:25]
* [-e] error rate               [default: 0.01]
* [-G] genome size              [default: 0]
* [-t] number of threads        [default: 1]
* [-o] output prefix file name  [default: LightAssembler]
``` 

#### Inputs 
LightAssembler assembles multiple input files of sequencing reads given in **fasta/fastq format**. Also, LightAssembler can read directly the input files compressed with gzip **fasta.gz/fastq.gz**.
