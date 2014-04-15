# How to get the EMu software

The current stable release and a pre-compiled executable of EMu for Mac OS X (64bit), can be found [here](https://github.com/andrej-fischer/EMu/releases).

# Test EMu

To test `EMu`, type the following on the command line while in the main directory.

`$ ./build/EMu --mut 21_breast_cancers.mutations.txt --opp 21_breast_cancers.opportunity.txt --pre ./test`

# Compilation

To compile EMu yourself, you need an installation of the GNU scientific library ([GSL](http://www.gnu.org/software/gsl/), v14 or later). On a machine without admin rights, you might need to change the paths in the Makefile to point to your local installation of GSL, if they are not in `/usr/local/`.
The software uses openMP. If you do not wish to use openMP, remove the `-fopenmp` flag in the Makefile.

Finally, simply type `make` on the command line.

# EMu usage

## Command line arguments: 

### Required

* `--mut [file]` Specify a file of mutation counts (a Nsamples x Nchannels matrix)

* `--opp [file/human-genome/human-exome]`  Specify the mutational opportunity.

  This can be either (i) the path to a flat text file of mutational opportunities (a Nsamples x Nchannels matrix) or (ii) `human-genome` to use the human whole genome opportunity for all samples or (iii) `human-exome` to use the whole human exome (female) for all samples.

* `--pre [path:./out]` Set the prefix used for all output files (e.g. `./here/results`).

### Optional

* `--force [int]`   Force the program to use a specific number of processes.

* `--mcmc [int]`    Run a MCMC with this number of steps to probe the posterior probability 
       		distribution for the mutational signatures and find error estimates.

* `--freeze [int]`  Perform zero-temperature Simulated-Annealing after convergence of the EM alorithm.

* `--spectra [file]` Use fixed mutational spectra (a Nspectra x Nchannels matrix).

    No EM will be performed. Only the activities per sample will be inferred and the mutations assigned. Useful for localizing processes in the genome.

* `--weights [file]` Supply (global) process activities to be used as an informed (local) activity prior. 
    This needs to be a (M x Nspectra) matrix, where Nsamples in `--mut` and `--opp` needs to be an integer multiple of M.
                  
## EMu output files:

* `^[pre]_[Nsp]_ml_spectra.txt`  The spectra found in the data using EM (Nspectra x Nchannels matrix)
* `^[pre]_[Nsp]_map_activities.txt`  The activities found in the data using EM (Nsamples x Nspectra matrix)
* `^[pre]_[Nsp]_assigned.txt`  The mutations assigned to each process (Nsamples x Nspectra matrix).
* `^[pre]_bic.txt`  The BIC values for the number of spectra tried.

### If MCMC was called:

* `^[pre]_[Nsp]_mcmc_spectra.txt`  The posterior mean spectra found in the data using MCMC (Nspectra x Nchannels matrix)
* `^[pre]_[Nsp]_mcmc_activities.txt`  The posterior mean activities found in the data using MCMC (Nsamples x Nspectra matrix)
* `^[pre]_[Nsp]_mcmc_err.txt`  The posterior std.dev. for the spectra using MCMC (Nspectra x Nchannels matrix)


# EMu-prepare usage

`EMu-prepare` is a program to create the input files for `EMu`.

## Command line arguments for `EMu-prepare`:

* `--mut [file]`  A flat text file with the mutations to be analysed. 

    Each line describes one mutation (please see note below). Expected format:
	 
	sample chromosome coordinate mutation
	
	sample: identifier for each sample (no white space)
	chomosome: integer (rename X=23,Y=24,mt=25 etc.)
	coordinate: one-based integer chomosome coordinate
	mutation: format A>T

* `--chr [dir]`	A directory of human chromosome fasta files.

    Expected file name format: `chr1.fa`. Rename file names for chr X, Y, mt etc., e.g. chrX.fa -> chr23.fa. You can download the latest version of the human reference genome [here](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/).

* `--cnv [file]`  A file with all the copy number information. 

    Each line is a non-standard copy number region. Format:

        sample chromosome start stop multiplier
	
    sample: identifier for each sample (no white space)
    chomosome: integer (rename X=23,Y=24,mt=25 etc.)
    start: chromosome start coordinate of cnv region
    stop: chromosome stop coordinate of cnv region (if -1, then extends to the end of the chromosome)
    multiplier: integer (in this region, this multiplier is used to integrate the opportunity)

    Note: the default multiplier is 2. This can be changed with --default [int]. If a sample has no copy number changes, still include at least one dummy line for each sample under consideration.

* `--pre [string]` A path for the bin-wise output files. 

    Since there will be one file for each sample and each chr, it is a good idea to send them to a separate directory.

* `--bin [int]` The size of the non-overlapping windows for which to get mutational/opportunity data.

* `--regions [file]`   A file with coordinates of sequenced genomic regions. 

    One region per line. Expected format:
        chromosome start stop

## EMu-prepare output files

Assuming `EMu-prepare` was called with `--cnv cnv.txt --mut mutations.txt`:

* `mutations.txt.96` The same as `mutations.txt`, with the mutation channel appended at the end of each line.

* `mutations.txt.mut.matrix`  A matrix of mutation counts with no. samples rows and 96 columns. Suitable for EMu.

* `mutations.txt.mut.samples`  The samples corresponding to each row in above file.

* `cnv.txt.opp.matrix`  A matrix of opportunity counts with no. samples rows and 96 columns. Suitable for EMu.

* `cnv.txt.opp.sample`  The samples corresponding to each row in above file. Check that this is the same order as in `mutations.txt.mut.samples`.

## NOTE:

In order to translate mutations to the 96 channels, `EMu-prepare` reads the bases 5' and 3' to the one given in a line of `mutations.txt` from the hard disk. It is very useful to sort the mutations file by chromosome and coordinate (otherwise the most time will be spent moving between physical locations in the hard disk). On UNIX, this can be achieved with:

`sort -k2n,2 -k3n,3 mutations.txt > mutations.sorted.txt`

# KNOWN BUGS

There is a known bug when openMP is compiled with the Mac OS compiler gcc version 4.2.1, which leads to random `abort trap:6` crashes. If possible, compile with latest gcc version. Alternatively, you can remove the `-fopenmp` flag from the Makefile or set the number of threads manually to one via:

`export OMP_NUM_THREADS=1; ./EMu --mut 21_breast_cancers.mutations --opp 21_breast_cancers.opportunity --pre ./target/test`
