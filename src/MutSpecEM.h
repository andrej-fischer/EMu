//MutSpecEM.h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// GSL headers...
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sf_gamma.h"

// parameter to the function to be minimized

struct Fpar{
  gsl_vector * cts;
  gsl_vector * ps_cts;
  gsl_matrix * M;
  int with_ps;
};


class MutSpecEM {
 private:
  gsl_matrix * spectra_start;
  gsl_matrix * counts;
  gsl_matrix * opp;
  gsl_matrix * pseudo_counts, * pseudo_counts_E, * pseudo_counts_M;
  gsl_matrix * overlap;
  gsl_matrix * weights_start;
  gsl_matrix * dist;
  gsl_vector * fexp;
  gsl_vector * tot_counts;
  gsl_vector * aff_chn;
  gsl_vector * done_with_ps;
  gsl_vector * factorials;
  gsl_vector * ps_const;
  int Nsp;// no. spectra
  int Nch;// no. channels
  int Nsa;// no. samples
  int Nopp;//no. opp spectra
  int no_vars;
  double spectra_min, weights_min;
  double total_opp;
  void set_weights_start();
  void set_dist();
  void set_pseudo_counts();
  void set_overlap();
  Fpar myFpar;
  
  void chop_weights();
  double get_mean_dev(); 
  double get_non_trivial_min(gsl_matrix * M);
  double new_dev;

 public:
  // Constructor
  MutSpecEM(  int Nspectra, gsl_matrix * Mut, gsl_matrix * Opp);
  //Destructor
  ~MutSpecEM();
  gsl_matrix * spectra;
  gsl_matrix * spectra_mean;
  gsl_matrix * spectra_err; 
  gsl_matrix * assigned_counts;  
  gsl_matrix * weights, * prior_weights;
  void pick_random_start();   
  int  get_weights_guess();
  int  get_spectra_guess();
  void get_all_map_weights( int loud, int with_ps);
  void get_all_map_spectra(int with_ps);
  double run_em(int loud, int max_em_steps, double& llh, double& bic);
  void get_error_estimates(int with_ps);
  void assign_mutations();
  //double random_search(int rs_steps, double llh);
  void run_mcmc(int mcmc_steps, const char * pre);
  void freeze_out(int steps);
  double get_llhood();
  double get_app_llhood(gsl_matrix * spectra_trial);
  void get_naive_ps_cts( gsl_matrix * Pseudo, gsl_matrix * Spectra, gsl_vector * Opp);
};



double get_map( gsl_vector * mix_sp, Fpar * myFpar);

double F_map(const gsl_vector * x, void * p);
void dF_map(const gsl_vector * x, void * p, gsl_vector * grad);
void FdF_map(const gsl_vector * x, void * p, double * val, gsl_vector * grad);

double get_pos_min(gsl_matrix * M);
double get_pos_min(gsl_vector * M);
void print_matrix(const gsl_matrix * M);
void print_vector(const gsl_vector * x);

double angle_to_simplex( gsl_vector * angle, gsl_vector * simplex);
void simplex_to_angle( gsl_vector * simplex, gsl_vector * angle);
void random_step_angle_uniform( gsl_vector * angle, double size);
void make_step_on_simplex( gsl_vector * simplex, double scale);
