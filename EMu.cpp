/* 

******************************************************************************

Copyright (c) 09-16-13  Genome Research Ltd.

Author: Andrej Fischer (af7[at]sanger.ac.uk)

This file is part of EMu v1.4

EMu is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation; either version 3 
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  
If not, see <http://www.gnu.org/licenses/>.


******************************************************************************

EMu.cpp

A program to infere mutational process signatures via the Expectation Maximization algorithm.

Command line arguments: 

Required:
--mut [file]    - The path to the flat text file of mutation counts (a Nsamples x Nchannels matrix)
--opp [file]    - The path to the flat text file of mutational opportunities (a Nsamples x Nchannels matrix)
--pre [path]    - The string to prefix the output files with (e.g. ./here/results)
Optional:
--force [int]   - Forces the program to use a specific number of processes for the fine search.

--mcmc [int]    - Performs a MCMC sampling with the desired number of steps to probe the posterior probability distribution
                  for the mutational signatures and fidn thus error estimates.
--freeze [int]  - Performs zero-temperature Simulated-Annealing after convergence of the EM alorithm.

--spectra [file] - Supply a matrix of mutational spectra (Nspectra x Nchannels matrix) to be used for the inference of activities.
                   No EM will be performed. Only activities per sample will be inferred and the mutations assigned.
--weights [file] - Supply a matrix of process activities to be used as an informed prior for the inference of activities per sample.
                   Needs to be a (M x Nspectra) matrix, where Nsamples in --mut and --opp needs to be a multiple of M.
                  
Output files:

^[pre]_[Nsp]_ml_spectra.txt      - The spectra found in the data using EM (Nspectra x Nchannels matrix)
^[pre]_[Nsp]_map_activities.txt  - The activities found in the data using EM (Nsamples x Nspectra matrix)
^[pre]_[Nsp]_assigned.txt        - The mutations assigned to each process (Nsamples x Nspectra matrix).
^[pre]_bic.txt                   - The BIC values for the number of spectra tried.

If MCMC was called:
^[pre]_[Nsp]_mcmc_spectra.txt     - The posterior mean spectra found in the data using MCMC (Nspectra x Nchannels matrix)
^[pre]_[Nsp]_mcmc_activities.txt  - The posterior mean activities found in the data using MCMC (Nsamples x Nspectra matrix)
^[pre]_[Nsp]_mcmc_err.txt         - The posterior std.dev. for the spectra using MCMC (Nspectra x Nchannels matrix)

*/

// General headers...
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>


// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"


// Own headers...
#include "MutSpecEM.h"

using namespace std;

// command line options...
struct cmdl_opts{
  const char * opp_fn;
  const char * mut_fn;
  const char * gspec_fn; 
  const char * gwts_fn;
  const char * prefix;
  int fNsp;
  int with_gspec, with_opp, with_gwts, mcmc, pseudo,freeze;
};
void get_opts( int argc, const char ** argv, cmdl_opts& opts);

// get the dimensions of the input data
void get_dims( const char * data_file, int * Ng, int * Ns);
// find the supported number of spectra in the data...
int find_no_spectra(
		    gsl_matrix * Mut, 
		    gsl_matrix * Opp, 
		    cmdl_opts& opts,
		    gsl_matrix ** nmax_spectra, 
		    gsl_matrix ** nmax_weights
		    );
// read the data into a gsl matrix, or print from...
void get_data( gsl_matrix * data, const char * data_file);
void print_data( gsl_matrix * data, const char * prefix, int Nsp, const char * suffix);

// ***MAIN***
int main(int argc, const char * argv[]){
  //
  int Nch = 1;// The number of mutation channels
  int Nsa = 1;// The number of samples
  int Nsp = 2;// The number of mutation spectra 
   // Preparing command line arguments...
  cmdl_opts opts;
  // Read options from the command line...
  get_opts( argc, argv, opts);
  // Get the input data dimensions...
  int Nch2=0,Nopp=0,Npat=0,Nbins=0;
  get_dims( opts.mut_fn, &Nsa, &Nch);
  // If opportunity was given explicitly...
  if ( opts.with_opp == 1){
    get_dims( opts.opp_fn, &Nopp, &Nch2);
    if (Nch != Nch2){
      cout<<"Dimensions "<<Nch<<" and "<<Nch2<<" are incompatible!"<<endl;
      exit(1);
    }
    if ( Nsa != Nopp && Nopp != 1){
      cout << "The no. samples " << Nsa << " and no. of opportunities "<<Nopp<< " are incompatible."<<endl;
      exit(1);
    }
  }  
  gsl_matrix * Mut = gsl_matrix_alloc( Nsa, Nch);
  gsl_matrix * Opp = gsl_matrix_alloc( Nsa, Nch);
  // Get the mutation data...
  get_data( Mut, opts.mut_fn);
  // Get the opportunity data, if given...
  if ( opts.with_opp == 1){
    if (Nopp == Nsa){
      get_data( Opp, opts.opp_fn);
    }
    // if only a single opportunity spectrum is given, apply to all..
    else if (Nopp == 1){
      gsl_matrix * opp = gsl_matrix_alloc( Nopp, Nch);
      get_data( opp, opts.opp_fn);
      for (int m=0; m<Nsa; m++){
	for (int j=0; j<Nch; j++){
	  gsl_matrix_set( Opp, m, j, gsl_matrix_get( opp, 0, j));
	}
      }
      gsl_matrix_free(opp);
    }
  }
  else{
    // Human Genome opportunity, if nothing was given...
    if (Nch==96){
      double HumOpp[96] = {
	1.14e+08,6.60e+07,1.43e+07,9.12e+07,//C>A @ AC[ACGT]
	1.05e+08,7.46e+07,1.57e+07,1.01e+08,//C>A @ CC[ACGT]
	8.17e+07,6.76e+07,1.35e+07,7.93e+07,//C>A @ GC[ACGT]
	1.11e+08,8.75e+07,1.25e+07,1.25e+08,//C>A @ TC[ACGT]
	1.14e+08,6.60e+07,1.43e+07,9.12e+07,//C>G @ AC[ACGT]
	1.05e+08,7.46e+07,1.57e+07,1.01e+08,//C>G @ CC[ACGT]
	8.17e+07,6.76e+07,1.35e+07,7.93e+07,//C>G @ GC[ACGT]
	1.11e+08,8.75e+07,1.25e+07,1.25e+08,//C>G @ TC[ACGT]
	1.14e+08,6.60e+07,1.43e+07,9.12e+07,//C>T @ AC[ACGT]
	1.05e+08,7.46e+07,1.57e+07,1.01e+08,//C>T @ CC[ACGT]
	8.17e+07,6.76e+07,1.35e+07,7.93e+07,//C>T @ GC[ACGT]
	1.11e+08,8.75e+07,1.25e+07,1.25e+08,//C>T @ TC[ACGT]
	1.17e+08,7.57e+07,1.04e+08,1.41e+08,//T>A @ AC[ACGT]
	7.31e+07,9.55e+07,1.15e+08,1.13e+08,//T>A @ CC[ACGT]
	6.43e+07,5.36e+07,8.52e+07,8.27e+07,//T>A @ GC[ACGT]
	1.18e+08,1.12e+08,1.07e+08,2.18e+08,//T>A @ TC[ACGT]
	1.17e+08,7.57e+07,1.04e+08,1.41e+08,//T>C @ AC[ACGT]
	7.31e+07,9.55e+07,1.15e+08,1.13e+08,//T>C @ CC[ACGT]
	6.43e+07,5.36e+07,8.52e+07,8.27e+07,//T>C @ GC[ACGT]
	1.18e+08,1.12e+08,1.07e+08,2.18e+08,//T>C @ TC[ACGT]
	1.17e+08,7.57e+07,1.04e+08,1.41e+08,//T>G @ AC[ACGT]
	7.31e+07,9.55e+07,1.15e+08,1.13e+08,//T>G @ AC[ACGT]
	6.43e+07,5.36e+07,8.52e+07,8.27e+07,//T>G @ AG[ACGT]
	1.18e+08,1.12e+08,1.07e+08,2.18e+08 //T>G @ AT[ACGT]
      };
      for (int m=0; m<Nsa; m++){
	for (int j=0; j<Nch; j++){
	  gsl_matrix_set( Opp, m, j, HumOpp[j]);
	}
      }
    }
    else{
      gsl_matrix_set_all( Opp, 1.0);
    }
  }
  // *** SCENARIO 1 ***
  // Read in externally supplied spectra and global weights to serve as prior...
  gsl_matrix * global_spectra = NULL;
  gsl_matrix * global_weights = NULL;  
  if ( opts.with_gwts == 1 && opts.with_gspec == 0){
    cout<<"Global weights can only be submitted WITH global spectra!\n";
    exit(1);
  }
  // if global weights are given, check global spectra...
  if ( opts.with_gspec == 1){
    if (opts.fNsp>0){
      cout<<"You cannot force a number of states with given spectra file!"<<endl;
      exit(1);
    }
    get_dims( opts.gspec_fn, &Nsp, &Nch2);
    if (Nch2 != Nch){
      printf("The provided spectra have %i channels instead of %i!\n",Nch2,Nch);
      exit(1);
    }
    global_spectra = gsl_matrix_alloc(Nsp,Nch);
    get_data( global_spectra, opts.gspec_fn);
    // Read in external global weights per sample...
    if ( opts.with_gwts == 1){
      get_dims( opts.gwts_fn, &Npat, &Nsp);
      if ( (Nsa % Npat) != 0 ){
	printf("The number of samples (%i) is not a multiple of the number of patients (%i).\n",Nsa,Npat);
	exit(1);
      }
      else{
	Nbins = (int) Nsa / Npat;
	printf("Found %i samples and %i global groups, i.e. %i bins each.\n",Nsa,Npat,Nbins);
      }
      global_weights = gsl_matrix_alloc(Npat,Nsp);
      get_data( global_weights, opts.gwts_fn);
    }    
  }
  // if global spectra and weights are supplied, we just need to find the (local) weights...
  if ( opts.with_gspec == 1){
    // Initiate model class...
    MutSpecEM myEM( Nsp, Mut, Opp);
    gsl_matrix_memcpy( myEM.spectra, global_spectra);
    if ( opts.with_gwts == 1){
      for (int m=0; m<Nsa; m++){
	int w = (int) m / Nbins;
	for (int a=0; a<Nsp; a++){
	  // set the activity weights to serve as priors...
	  gsl_matrix_set( myEM.prior_weights, m, a, gsl_matrix_get( global_weights, w, a));
	}
      }
    }
    myEM.get_weights_guess();
    int loud=1;
    // infer the weights for the supplied samples...
    myEM.get_all_map_weights( loud, opts.pseudo);
    //double llh = myEM.get_llhood();
    //printf("\nLLH = %e\n", llh);
    myEM.assign_mutations();
    print_data( myEM.weights, opts.prefix, Nsp, "activities");
    print_data( myEM.assigned_counts, opts.prefix, Nsp, "assigned");
    gsl_matrix_free( global_spectra);
    if ( opts.with_gwts == 1) gsl_matrix_free( global_weights);
    cout<<endl;
  }
  // *** SCENARIO 2 ***
  // If nothing is supplied, we need to infer the spectra...
  else{
    int nmax=1;
    gsl_matrix * nmax_spectra = NULL;
    gsl_matrix * nmax_weights = NULL;
    //struct timeval t1, t2;
    //gettimeofday( &t1, NULL);
    // if no. spectra is forced, use that...
    if (opts.fNsp>0){
      printf("Will use %i spectra, as requested.\n", opts.fNsp);
    }
    else{
      // *** QUICK SCREEN *** 
      // (over Nsp to find the one with max bic...)
      nmax = find_no_spectra(Mut,Opp,opts, &nmax_spectra, &nmax_weights);
    }
    // set the no. spectra to the one with maximal BIC or the forced one... 
    Nsp = (opts.fNsp>0) ? opts.fNsp : nmax;
    gsl_matrix * best_spectra  = gsl_matrix_alloc(Nsp,Nch);
    gsl_matrix * best_weights  = gsl_matrix_alloc(Nsa,Nsp);
    MutSpecEM myEM( Nsp, Mut, Opp);
    // if forced Nsp, try 10 random starting positions and pick the best...
    if (opts.fNsp > 0){
      double l,b,maxb;
      int trials=10;
      for (int t=0; t<trials; t++){
	printf("\rRandom trial %i of %i...",t+1,trials);
	cout<<flush;
	myEM.pick_random_start();
	myEM.run_em( 0, 0, l, b);
	if (t==0 || b>maxb){
	  gsl_matrix_memcpy( best_spectra, myEM.spectra);
	  gsl_matrix_memcpy( best_weights, myEM.weights);
	  maxb = b;
	}
      }
      cout<<endl;
    }
    else{
      gsl_matrix_memcpy( best_spectra, nmax_spectra);
      gsl_matrix_memcpy( best_weights, nmax_weights);
      gsl_matrix_free(nmax_spectra);
      gsl_matrix_free(nmax_weights);
    }
    // re-insert the winners of the random trials above to be followed through deep EM...
    gsl_matrix_memcpy( myEM.spectra, best_spectra);
    gsl_matrix_memcpy( myEM.weights, best_weights);
    // print all data before EM...
    print_data( myEM.spectra, opts.prefix, Nsp, "ml_spectra");
    print_data( myEM.weights, opts.prefix, Nsp, "map_activities");
    // *** RUN DEEP EM HERE ***
    int loud = 1, num_em = 1.0e4;
    double llh,bic;
    myEM.run_em( loud, num_em, llh, bic);
    if (opts.freeze>0) myEM.freeze_out(opts.freeze);
    // re-print data after EM (and freeze-out)...
    print_data( myEM.spectra, opts.prefix, Nsp, "ml_spectra");
    print_data( myEM.weights, opts.prefix, Nsp, "map_activities");
    myEM.get_error_estimates(1);
    print_data( myEM.spectra_err, opts.prefix, Nsp, "spectra_ErrEst");
    // get the mutation assignments...
    loud = 1;
    int with_ps = 1;
    myEM.get_all_map_weights( loud, with_ps);
    myEM.assign_mutations();
    print_data( myEM.assigned_counts, opts.prefix, Nsp, "assigned");
    // *** DO MCMC HERE ***
    if (opts.mcmc>0){
      myEM.run_mcmc(opts.mcmc,opts.prefix);
      print_data( myEM.spectra_mean, opts.prefix, Nsp, "mcmc_spectra");
      print_data( myEM.spectra_err, opts.prefix, Nsp, "mcmc_err");
      print_data( myEM.weights, opts.prefix, Nsp, "mcmc_activities");
    }    
    /*
    // print the EM search execution time
      gettimeofday(&t2,NULL);
      // Calculate time it took
      //double accum = ( t2.tv_sec - t1.tv_sec ) + ( t2.tv_usec - t1.tv_usec ) / 1.0e6;
      // printf( "%lf\n", accum );
      char buff[128];
      sprintf( buff,"%s_%i_time.txt", opts.prefix, Nsp);
      string time_fn = buff;
      FILE * outfile_pt = fopen( time_fn.c_str(), "w");
      //fprintf( outfile_pt, "%.3f\n", (double) (t2 - t1)/CLOCKS_PER_SEC);
      fprintf( outfile_pt, "%.6f\n",accum);
      fclose(outfile_pt);
    */
    // cleanup...   
    gsl_matrix_free( best_weights);
    gsl_matrix_free( best_spectra);
  }
  gsl_matrix_free( Opp);
  gsl_matrix_free( Mut);
  //
  cout<<endl;
  return(0);
}



//*** FUNCTIONS ***

// Get the size of the data matrices...
void get_dims(const char * infile_name, int * dim1, int * dim2){
  ifstream inf_str;
  inf_str.open(infile_name);
  string line;
  double val=0;
  int line_ct=0, col_ct=0, old_col_ct=0;
  while ( getline(inf_str, line) ){
    line_ct++;
    stringstream line_ss;
    line_ss.str(line);
    while( line_ss >> val ){
      col_ct++;
    }
    if(old_col_ct==0){
      old_col_ct = col_ct;
    }
    else if(col_ct != old_col_ct){
      printf("Found a different number of columns (%i) in line %i!\n", col_ct, line_ct);
      exit(1);
    }
    col_ct = 0;    
  }
  *dim1 = line_ct;
  *dim2 = old_col_ct;
  //printf("Found %i rows and %i columns.\n", line_ct, old_col_ct);
  inf_str.close();
}

// Read data from a flat text file into gsl-matrix...
void get_data(gsl_matrix * data, const char * infile_name){
  FILE * data_fp = fopen( infile_name, "r");
  int err = 1;
  err = gsl_matrix_fscanf(data_fp, data);
  fclose(data_fp);
  if ( err != 0){
    printf("reading not successful!\n");
    exit(1);
  }
}
 

void print_data( gsl_matrix * data, const char * prefix, int Nsp, const char * suffix){
  char buff[128];
  string  out_fn;
  sprintf( buff,"%s_%i_%s.txt", prefix, Nsp, suffix);
  out_fn = buff;
  FILE * outfile_pt = fopen( out_fn.c_str(), "w");
  for (int i = 0; i< (int) data->size1; i++){
    for (int j = 0; j< (int) data->size2; j++){
      fprintf( outfile_pt, "%.3e ", gsl_matrix_get( data, i, j));
    }
    fprintf( outfile_pt, "\n");
  }
  fclose(outfile_pt);    
}



// get command line arguments...
void get_opts( int argc, const char ** argv, cmdl_opts& opts){
  // prepare...
  opts.mut_fn   = NULL;
  opts.opp_fn   = NULL;
  opts.gspec_fn = NULL;
  opts.gwts_fn  = NULL;
  opts.prefix = "out";
  opts.with_gspec = 0;
  opts.with_gwts  = 0;
  opts.with_opp   = 0;
  opts.fNsp   = 0; 
  opts.mcmc   = 0;
  opts.freeze   = 0;
  opts.pseudo = 0;
  // get the command line arguments...
  int opt_ind = 1;
  while( opt_ind < argc && (argv[opt_ind][0] == '-') ){ 
    string opt_switch = argv[opt_ind];
    if(opt_switch == "--force"){
      opt_ind++;
      opts.fNsp = atoi( argv[opt_ind] );
    }
    else if(opt_switch == "--mcmc"){
      opt_ind++;
      opts.mcmc = (int) atof( argv[opt_ind] );
    }
    else if(opt_switch == "--freeze"){
      opt_ind++;
      opts.freeze = (int) atof( argv[opt_ind] );
    }
    else if(opt_switch == "--pseudo"){
      opts.pseudo = 1;
    }
    else if (opt_switch == "--mut"){
      opt_ind++;
      opts.mut_fn = argv[opt_ind];
    }
    else if (opt_switch == "--opp"){
      opt_ind++;
      opts.opp_fn = argv[opt_ind];
      opts.with_opp = 1;
    }
    else if (opt_switch == "--spectra"){
      opt_ind++;
      opts.gspec_fn = argv[opt_ind];
      opts.with_gspec = 1;
    }
    else if (opt_switch == "--weights"){
      opt_ind++;
      opts.gwts_fn = argv[opt_ind];
      opts.with_gwts = 1;
      opts.pseudo = 1;
    }
    else if (opt_switch == "--pre"){
      opt_ind++;
      opts.prefix = argv[opt_ind];
    }
    else if (opt_switch == "--version"){
      opt_ind++;
      cout<<"This is version 1.2 of EMu, an EM-algorithm for the inference of mutational signatures.\n";
      exit(0);
    }
    else{
      cout<<"Usage: ./EMu --mut 21_breast_cancers.mutations --opp 21_breast_cancers.opportunity --pre ./target/test"<<endl;
      exit(1);
    }
    opt_ind++;
  }
}

// Determining the number of spectra in the data via BIC...
int find_no_spectra( 
		     gsl_matrix * Mut, 
		     gsl_matrix * Opp, 
		     cmdl_opts& opts,
		     gsl_matrix ** nmax_spectra, 
		     gsl_matrix ** nmax_weights 
		      ){
  int Nsa = (int) Mut->size1;
  int Nch = (int) Mut->size2;
  int nmax=0,was_close=0;
  // remember the best solutions in each case...
  gsl_matrix **  best_trial_spectra = new gsl_matrix * [Nsa];
  gsl_matrix **  best_trial_weights = new gsl_matrix * [Nsa];
  int * allocated = new int [Nsa];
  double * bics = new double [Nsa];
  for (int m=0; m<Nsa; m++){
    allocated[m]=0;
  }
  double l,b,bic,max_bic,llh;
  for (int ntrial=1; ntrial<Nsa; ntrial++){
    (ntrial==1) 
      ? printf("Trying to find %2i spectrum...", ntrial) 
      : printf("Trying to find %2i spectra....", ntrial);
    // Initiate Model class...
    MutSpecEM trialEM( ntrial, Mut, Opp);
    best_trial_spectra[ntrial] = gsl_matrix_alloc(ntrial,Nch);
    best_trial_weights[ntrial] = gsl_matrix_alloc(Nsa,ntrial);
    allocated[ntrial] = 1;
    for (int t=0; t<10;t++){
      int loud=0, num_em=0;
      trialEM.pick_random_start();
      trialEM.run_em( loud, num_em, l, b);
      if (t==0 || b > bic){
	bic = b;
	llh = l;
	gsl_matrix_memcpy(best_trial_spectra[ntrial],trialEM.spectra);
	gsl_matrix_memcpy(best_trial_weights[ntrial],trialEM.weights);
      }
    }
    printf("LLH = %.6e BIC = %.6e\n", llh, bic);
    bics[ntrial] = bic;
    // test for stopping criterion...
    if (ntrial>1 && fabs(bic-max_bic)<20.0){// close calls...
      was_close=1;
      cout<<"Close call!\n";
      cout<<"\nRe-fining EM for n = "<<ntrial-1<<endl;
      // do deep EM for last...
      MutSpecEM redoEM( ntrial-1, Mut, Opp);
      gsl_matrix_memcpy(redoEM.spectra, best_trial_spectra[ntrial-1]);
      gsl_matrix_memcpy(redoEM.weights, best_trial_weights[ntrial-1]);
      redoEM.run_em( 1, 1000, l, max_bic);
      cout<<"\nRe-fining EM for n = "<<ntrial<<endl;
      // do deep EM for this...
      gsl_matrix_memcpy(trialEM.spectra, best_trial_spectra[ntrial]);
      gsl_matrix_memcpy(trialEM.weights, best_trial_weights[ntrial]);
      trialEM.run_em( 1, 1000, l, bic);	
    }
    if (ntrial>1  && (bic <= max_bic || was_close==1)){
      break;
    }
    else if (bic > max_bic || ntrial == 1){
      max_bic = bic;
      nmax = ntrial;
    }
  }
  printf("\nAccording to BIC, there are %i spectra in the data.\n\n",nmax);
  // print bics to file...
  char buff[128];
  sprintf( buff,"%s_bic.txt", opts.prefix);
  string bic_fn = buff;
  FILE * outfile_pt = fopen( bic_fn.c_str(), "w");
  for (int n=1; n<=nmax+1; n++){
    fprintf( outfile_pt, "%i %e\n", n, bics[n]);
  }
  fclose(outfile_pt);
  //export the best solutions for the supported no. processes...
  *nmax_spectra = gsl_matrix_alloc(nmax,Nch);
  *nmax_weights = gsl_matrix_alloc(Nsa,nmax);
  gsl_matrix_memcpy(*nmax_spectra,best_trial_spectra[nmax]);
  gsl_matrix_memcpy(*nmax_weights,best_trial_weights[nmax]);
  // then free the registers...
  for (int m=1; m<Nsa; m++){
    if (allocated[m] == 1){
      gsl_matrix_free(best_trial_spectra[m]);
      gsl_matrix_free(best_trial_weights[m]);
    }
  }
  delete [] allocated;
  delete [] bics;
  delete [] best_trial_spectra;
  delete [] best_trial_weights;
  return(nmax);
}
