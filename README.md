## LightAssembler
Lightweight resources assembly algorithm for high-throughput sequencing reads.
#### System requirements 
64-bit machine, g++ compiler or gcc in general, [pthreads](http://en.wikipedia.org/wiki/POSIX_Threads),and [zlib](http://en.wikipedia.org/wiki/Zlib) libraries.

#### Installation 
1. Clone the [GitHub repo](https://github.com/SaraEl-Metwally/LightAssembler), e.g. with `git clone https://github.com/SaraEl-Metwally/LightAssembler.git`
2. Run `make` in the repo directory for **k <= 31**  or `make k=kmersize` for **k > 31**, e.g. `make k=49`. 

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

#### Notes
- If the gap size parameter is missing, LightAssembler invokes its parameters extrapolation module to compute the starting gap based on the sequencing coverage and the error rate of the dataset.
- The maximum read length for this version is ``` 1024```.
- The maximum supported read files for this version ```100```.

#### Input read files 
LightAssembler assembles multiple input files of the sequencing reads given in ***fasta/fastq*** format. Also, LightAssembler can read directly the input files compressed with gzip ***fasta.gz/fastq.gz***.

#### Outputs
The output of LightAssembler is the set of assembled contigs in ***fasta*** format, in the file:

``` [output prefix].contigs.fasta ``` 

LightAssembler also reports the following on the screen:

- Number of resulted contigs.
- Maximum contig length.
- Total Assembly size.
- Total genome coverage.
- Total Assembly time as well as the total time for each step.

Also, by using the ```--verbose``` option, LightAssembler reports additional details for each step such as the number of kmers, the false positive rate of Bloom filter and the number of branching kmers in the dataset, average read length and the average sequencing coverage.

#### Example 1
``` ./LightAssembler -k 31 -g 15 -e 0.01 -G 4686137 -o ecoli_contigs -t 3 ecoli_reads_1.fq ecoli_reads_2.fq --verbose ```

```
--- Uniform kmers sampling. 

--- h(0):m(0):s(5) elapsed time.
--- total number of kmers in BloomA = 7791111
--- BloomA false positive rate = 0.00193375
--- average read length = 101
--- average sequencing coverage = 35
--- probability of an incorrect kmer appears in the sample : 0.0249524

--- Trusted/untrusted kmers filtering. 

--- h(0):m(0):s(24) elapsed time.
--- total number of kmers in BloomB = 4548112
--- BloomB false positive rate = 7.7715e-05

--- Branching-kmers computation. 

--- h(0):m(0):s(5) elapsed time.
--- number of branching kmers = 54644

--- Graph traversal. 

--- h(0):m(0):s(16) elapsed time.
--- number of contigs     = 731
--- maximum contig length = 120924
--- assembly size         = 4473869
--- genome coverage       = 95.4703%

--- The assembly session is finished. 

--- h(0):m(0):s(31) elapsed time.```

#### Example 2 (missing -g )
