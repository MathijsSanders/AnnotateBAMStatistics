# AnnotateBAMStatistics

This repository contains the code of AnnotateBAMStatistics. AnnotateBAMStatistics is a multi-threaded algorithm for annotating ANNOVAR output with additional useful fragment count statistics for one ore more specified BAM files. This includes fragment coverage, number of mutant fragments, fragment-based variant allele frequency and mutant read bias for high quality and all fragments separately.

## How do I run it?

AnnotateBAMStatistics takes an ANNOVAR file and one or more BAM files (comma-separated list) as input and adds informative fragment-based statistics. These statistics are generally used for filtering the variant list or to determine clonality in cancer samples. AnnotateBAMStatistics writes the output to standard output.

### The easiest way

The application has been containerized into a Singularity container and is available from Singularity hub. In case singularity is installed, simply run:

```bash
singularity pull shub://MathijsSanders/AnnotateBAMStatisticsSingularity
```

The following parameters are available:

- -a/--annovarfile*: Input ANNOVAR file
- -b/--bamfiles*: One or more BAM files for which fragment statistics are calculated (comma-separated)
- -m/--min-alignment-score: The minimum alignment score threshold to determine whether a read is of high quality (Default: 40, based on BWA)
- -t/--threads: Number of threads (Default 1)

### The more difficult way

First clone the GitHub repository with git and run the following line:

```bash
make all
```
The binary AnnotateBAMStatistics is located at:

```bash
./dist/Release/GNU-Linux/AnnotateBAMStatistics
```

The following parameters are available:

- -a/--annovar-file*: Input ANNOVAR file
- -b/--bam-files*: One or more BAM files for which fragment statistics are calculated (comma-separated)
- -s/--min-base-score: Minimum base score for filtered statistics (Default=30)
- -S/--min-alignment-score: Minimum alignment score for filtered statistics (Default=40, based on BWA)
- -m/--min-match-length: Minimum alignment length for a read to be considered (Default=10)
- -t/--threads: Number of threads (Default 1)
- -d/--count-duplicates: If duplicates are used in the statistics (Default=false)
- -u/--count-supplementary: If supplementary alignments are used in the statistics (Default=false)
- -r/--pileup-regions: If individual regions from the annovar files are piled up. Otherwise scan the genome (Default=false)
- -v/--verbose: If specified be verbose (Default=false)
- -h/--help: Help information

*Dependencies*
- gcc 4.3+
- g++ 4.3+
