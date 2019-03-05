# README #

###Metagenomic Analysis###

Unique marker genes from MetaPhlAn2 are used to detect the presence/absence of viral, bacterial and archaea strains and subsequently, species counts are estimated using kallisto

The following sections contain tutorials on how to run the D1 and D2 pipelines as described in

Teo, Y. V. & Neretti, N. A comparative study of metagenomics analysis pipelines at the species level. bioRxiv (2016).

###Quick setup###

**Dependencies**

The following versions are tested for the pipeline:

samtools/1.3.1

kallisto/0.42.2

bwa/0.7.12

```
module load samtools/1.3.1
module load kallisto/0.42.2
module load bwa/0.7.12
```

###How to run?###

KallistoMarkerIdx, FullGenome and Species List files are available upon request at yee_voan_teo[AT]brown.edu

**For samples from human-associated habitat**

```
usage: BK.py [-h] [-o O] [-t T] [-infile1 INFILE1] [-infile2 INFILE2]
             [-outputprefix OUTPUTPREFIX]
             [-KallistoMarkerIdx KALLISTOMARKERIDX] [-FullGenome FULLGENOME]
             [-SpeciesList SPECIESLIST] [-bwaindex BWAINDEX]

Metagenomic Analysis of bacterial, viral and archael communities at the
species level from human-associated habitat

optional arguments:
  -h, --help            show this help message and exit
  -o O                  Output directory
  -t T                  Number of threads [1]
  -infile1 INFILE1      Path to gzip fastq 1
  -infile2 INFILE2      Path to gzip fastq 2
  -outputprefix OUTPUTPREFIX
                        Output prefix
  -KallistoMarkerIdx KALLISTOMARKERIDX
                        Kallisto markers index
  -FullGenome FULLGENOME
                        Full genome fasta
  -SpeciesList SPECIESLIST
                        Path to species list
  -bwaindex BWAINDEX    Path to BWA human index
```

**For samples from nonhuman-associated habitat**
```
usage: K.py [-h] [-o O] [-t T] [-infile1 INFILE1] [-infile2 INFILE2]
            [-outputprefix OUTPUTPREFIX]
            [-KallistoMarkerIdx KALLISTOMARKERIDX] [-FullGenome FULLGENOME]
            [-SpeciesList SPECIESLIST]

Metagenomic Analysis of bacterial, viral and archael communities at the
species level from NON-human associated habitat

optional arguments:
  -h, --help            show this help message and exit
  -o O                  Output directory
  -t T                  Number of threads [1]
  -infile1 INFILE1      Path to gzip fastq 1
  -infile2 INFILE2      Path to gzip fastq 2
  -outputprefix OUTPUTPREFIX
                        Output prefix
  -KallistoMarkerIdx KALLISTOMARKERIDX
                        Kallisto markers index
  -FullGenome FULLGENOME
                        Full genome fasta
  -SpeciesList SPECIESLIST
                        Path to species list
```

**Examples on how to run**

For samples that consist of human reads
```
/home/Metagenome/BK.py -o /home/Metagenome/test/ -outputprefix test -infile1 /home/Metagenome/test/test_1.fq.gz -infile2 /home/Metagenome/test/test_2.fq.gz \
-KallistoMarkerIdx /home/Metagenome/data/Marker_kallistoIdx -FullGenome /home/Metagenome/data/MicrobialCompleteGenome.fa \
-SpeciesList /home/Metagenome/data/Microbial_ID -bwaindex /data/hg19.fa
```

For samples that do not have human reads

```
/home/Metagenome/K.py -o /home/Metagenome/test/ -outputprefix test -infile1 /home/Metagenome/test/test_1.fq.gz -infile2 /home/Metagenome/test/test_2.fq.gz \
-KallistoMarkerIdx /home/Metagenome/data/Marker_kallistoIdx -FullGenome /home/Metagenome/data/MicrobialCompleteGenome.fa \
-SpeciesList /home/Metagenome/data/Microbial_ID 
```

###Output files###

Final species counts are in $OUTPUTPREFIX_FinalCount.txt file. Example output:
```
Shewanella_loihica  414.323
Rickettsia_canadensis   106.819
Human_herpesvirus_7 14.0
```
kallisto output at the strain level is in /$OUTPUTPREFIX/abundance.tsv