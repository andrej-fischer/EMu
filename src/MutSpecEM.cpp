//MutSpecEM.cpp

#define PI 3.1415926

// Own headers...
#include "MutSpecEM.h"

using namespace std;

// Constructor...
MutSpecEM::MutSpecEM( int Nspectra, gsl_matrix * Counts, gsl_matrix * Opp){
  // Read the data dimensions...
  Nsp = Nspectra;
  Nch = (int) Counts->size2;
  Nsa = (int) Counts->size1;
  // the number of parmeters of the model...
  no_vars = Nsp*(Nch-1) + Nsp*Nsa;
  // Allocate the matrices
  weights       = gsl_matrix_calloc(Nsa,Nsp);//
  weights_start = gsl_matrix_calloc(Nsa,Nsp);//
  prior_weights = gsl_matrix_calloc(Nsa,Nsp);//
  gsl_matrix_set_all( prior_weights, 1.0);// un-informative prior
  // data copy...
  counts = gsl_matrix_calloc( Nsa, Nch);
  assigned_counts = gsl_matrix_calloc( Nsa, Nsp);
  overlap = gsl_matrix_calloc( Nsa, Nsp);
  pseudo_counts   = gsl_matrix_calloc( Nsa, Nsp*Nch);
  pseudo_counts_E = gsl_matrix_calloc( Nsa, Nsp);
  pseudo_counts_M = gsl_matrix_calloc( Nch, Nsp);
  done_with_ps = gsl_vector_calloc(Nsa);
  opp = gsl_matrix_calloc(Nsa,Nch);// the opportunity
  if (gsl_matrix_isnonneg(Counts) == 1){
    gsl_matrix_memcpy( counts, Counts);
  }
  else{
    printf("Check mut input for negative entries!\n");
    exit(1);
  }
  if (gsl_matrix_isnonneg(Opp) == 1){
    gsl_matrix_memcpy( opp, Opp);
  }
  else{
    printf("Check opp input for negative entries!\n");
    exit(1);
  }
  // check and correct consistency: where no-opp, there no-counts!
  int screamed=0;
  for (int m=0;m<Nsa;m++){
    for (int j=0; j<Nch; j++){
      if ( gsl_matrix_get(opp,m,j) == 0 && gsl_matrix_get(counts,m,j) > 0){
	if (screamed==0){
	  cout<<"Found mutations where there should be none (opp=0). Set to zero!"<<endl;
	  screamed=1;
	}
	gsl_matrix_set(counts,m,j,0);
      }
    }
  }
  // get the total counts per sample...
  tot_counts = gsl_vector_calloc(Nsa);
  aff_chn = gsl_vector_calloc(Nsa);
  double tot;
  for (int m=0; m<Nsa; m++){
    tot=0;
    for (int j=0; j<Nch; j++){
      tot += gsl_matrix_get( counts, m, j);
      if (gsl_matrix_get( counts, m, j) > 0) (aff_chn->data[m])++;
    }
    gsl_vector_set(tot_counts,m,tot);
  }
  // set the factorial factors (stay const., go into log-likelihood)...
  factorials = gsl_vector_calloc(Nsa);
  ps_const   = gsl_vector_calloc(Nsa);
  for (int m=0; m<Nsa; m++){
    tot=0;
    for (int j=0; j<Nch; j++){
      tot += gsl_sf_lngamma( gsl_matrix_get( counts, m, j) + 1);
    }
    gsl_vector_set(factorials,m,tot);
  }
  // spectra matrices...
  spectra = gsl_matrix_calloc(Nsp,Nch);
  spectra_start = gsl_matrix_calloc(Nsp,Nch);
  spectra_mean = gsl_matrix_calloc(Nsp,Nch);
  spectra_err = gsl_matrix_calloc(Nsp,Nch);
  dist = gsl_matrix_alloc(Nsa,Nch);
  spectra_min = (double) 1.0e-6 / Nch;
  weights_min = 0.1;
  fexp = gsl_vector_calloc(Nsa);
  srand(time(NULL)); 
  //srand(1); 
  MutSpecEM::pick_random_start();
}



// Destructor...
MutSpecEM::~MutSpecEM(){
  gsl_matrix_free( spectra);
  gsl_matrix_free( spectra_start);
  gsl_matrix_free( spectra_mean);
  gsl_matrix_free( spectra_err);
  gsl_matrix_free( counts);
  gsl_matrix_free( pseudo_counts);
  gsl_matrix_free( pseudo_counts_E);
  gsl_matrix_free( pseudo_counts_M);
  gsl_matrix_free( assigned_counts);
  gsl_matrix_free( weights);
  gsl_matrix_free( prior_weights);
  gsl_matrix_free( weights_start);
  gsl_matrix_free( overlap);
  gsl_matrix_free( dist);  
  gsl_vector_free( fexp);
  gsl_matrix_free( opp);
  gsl_vector_free( tot_counts);
  gsl_vector_free( aff_chn);
  gsl_vector_free( factorials);
  gsl_vector_free( ps_const);
}



// set the entries of the emit_mean matrix to the observation means...
void MutSpecEM::pick_random_start(){
  gsl_vector_view  spec;
  for (int a=0; a<Nsp; a++){
    for (int j=0; j<Nch; j++){
      gsl_matrix_set( spectra, a, j, (double) rand() / RAND_MAX);
    }
  } 
  // normalize spectra
  for (int a=0; a<Nsp; a++){
    spec = gsl_matrix_row( spectra, a);
    gsl_vector_scale( &spec.vector, 1.0 / gsl_blas_dasum(&spec.vector) );
  }
  // set uniform weights...
  gsl_matrix_set_all( weights, 1.0e-3);
}


void  MutSpecEM::set_overlap(){
  double ol = 0;
  for (int m=0; m<Nsa; m++){
    for (int a=0; a<Nsp; a++){
      gsl_vector_view spectrum = gsl_matrix_row( spectra, a);
      gsl_vector_view oppV = gsl_matrix_row( opp, m);
      gsl_blas_ddot( &spectrum.vector, &oppV.vector, &ol);
      gsl_matrix_set( overlap, m, a, ol);
    }
  }
}

void MutSpecEM::set_pseudo_counts(){
  gsl_matrix_set_zero(pseudo_counts);
  gsl_matrix_set_zero(pseudo_counts_E);
  gsl_matrix_set_zero(pseudo_counts_M);
  double val, pw;
  gsl_vector * ps_cts = gsl_vector_alloc(Nsp*Nch);
  for (int m=0; m<Nsa; m++){
    gsl_vector_view oppV = gsl_matrix_row(opp,m);
    if (gsl_blas_dasum(&oppV.vector) == 0) continue;
    for (int a=0; a<Nsp; a++){
      pw = gsl_matrix_get( prior_weights, m, a);
      for (int j=0; j<Nch; j++){
	val  = gsl_matrix_get(spectra, a, j);
	val *= gsl_matrix_get(opp, m, j);
	val *= pw;
	gsl_vector_set(ps_cts, a*Nch + j, val);
      }
    }
    val = gsl_blas_dasum(ps_cts);
    if (val>0){
      gsl_vector_scale( ps_cts, (double) Nsp / val);
      gsl_matrix_set_row( pseudo_counts, m, ps_cts);
    }
    else{
      cout<<"ERROR in set_pseudo_counts(): val = "<<val<<endl;
      exit(1);
    }
  }
  gsl_vector_free(ps_cts);
  // the pseudo-counts for the E-step...
  for (int m=0; m<Nsa; m++){
    for (int a=0; a<Nsp; a++){
      val = 0;
      for (int j=0; j<Nch; j++){
	val += gsl_matrix_get( pseudo_counts, m, a*Nch + j);
      }
      gsl_matrix_set( pseudo_counts_E, m, a, val);
    }
  }
  // the pseudo-counts for the M-step...
  for (int j=0; j<Nch; j++){
    for (int a=0; a<Nsp; a++){
      val = 0;
      for (int m=0; m<Nsa; m++){
	val += gsl_matrix_get( pseudo_counts, m, a*Nch + j);
      }
      gsl_matrix_set( pseudo_counts_M, j, a, val);
    }
  }
  // the constant part of the log-likelihood...
  gsl_vector_set_zero(ps_const);
  double sum;
  for (int m=0; m<Nsa; m++){
    sum=0;
    gsl_vector_view ctsV = gsl_matrix_row( pseudo_counts, m);
    if ( gsl_blas_dasum(&ctsV.vector) == 0) continue;
    for (int a=0; a<Nsp; a++){
      for (int j=0; j<Nch; j++){
	val  = gsl_matrix_get( pseudo_counts, m, a*Nch + j);
	if (val>0) val *= log( gsl_matrix_get(spectra, a, j) * gsl_matrix_get(opp, m, j) );
	sum += val;
      }
    }
    gsl_vector_set( ps_const, m, sum);    
  }
}


void MutSpecEM::assign_mutations(){
  MutSpecEM::set_overlap();
  gsl_matrix_set_zero(assigned_counts);
  gsl_vector * ass_cts = gsl_vector_alloc(Nsp);
  for (int m=0; m<Nsa; m++){
    gsl_vector_view oppV = gsl_matrix_row( opp, m);
    if (gsl_blas_dasum(&oppV.vector) == 0 || gsl_vector_get(tot_counts,m) == 0){
      continue;
    }
    for (int a=0; a<Nsp; a++){
      gsl_vector_set( ass_cts, a, gsl_matrix_get( overlap, m, a) * gsl_matrix_get( weights, m, a));
    }
    gsl_vector_scale( ass_cts, gsl_vector_get(tot_counts,m) / gsl_blas_dasum(ass_cts));
    gsl_matrix_set_row( assigned_counts, m, ass_cts);
  }  
  gsl_vector_free(ass_cts);
}



// the EM routine for fixed starting conditions (in emit_start)...
double MutSpecEM::run_em(int loud, int max_em_steps, double& rllh, double& rbic){
  int em_steps = 0, max_app_steps = 1.0e3, is_sing=0;
  double test = 10, new_llh = 10, old_llh = 0, old_dev=0;
  // find the new mixtures starting from equal mix...
  int r=0;
  int em_loud=0, with_ps = 1;
  // approximative EM-Iteration...
  while ( r < max_app_steps ){
    // E-STEP
    is_sing = MutSpecEM::get_weights_guess();
    if (is_sing==1){
      //cout << "Singular matrix in weights_guess!\n";
      break;
    }
    // M-STEP
    is_sing = MutSpecEM::get_spectra_guess();
    if (is_sing==1) {
      //cout << "Singular matrix in spectra_guess!\n";
      break;
    }
    new_dev = MutSpecEM::get_mean_dev();
    test = (old_dev - new_dev) / new_dev;
    old_dev = new_dev;
    r++;
    if (loud ==1) printf("\rEM step %5i test %+.3e dev %.2e", r, test, new_dev);
    cout<<flush;
    if ( r > 1 && test < 1.0e-6 ) break;
  }
  if (loud ==1) cout<<endl;
  if (is_sing==1) max_em_steps = 1.0e3;
  if ( max_em_steps == 0){
    // E-STEP: First find the most likely activities 
    MutSpecEM::get_all_map_weights( em_loud, with_ps);
    // compute the data log-likelihood...
    new_llh = MutSpecEM::get_llhood();
  }
  else{
    // THIS IS THE FULL EM-ITERATION...
    while(  em_steps < max_em_steps){
      // E-STEP: First find the most likely activities 
      MutSpecEM::get_all_map_weights( em_loud,with_ps);
      // compute the data log-likelihood...
      new_llh = MutSpecEM::get_llhood();
      // compute test observables...
      new_dev = MutSpecEM::get_mean_dev();
      test =  (new_llh - old_llh) / fabs(new_llh);
      old_dev = new_dev;
      em_steps++;
      // report...
      if (loud == 1){
	printf("\rEM step %5i llh %.6e test %+.3e dev %.5e", em_steps, new_llh, test, new_dev);
	cout<<flush;
      }
      // stop condition...
      if ( em_steps > 1 && (test < 1.0e-7 || fabs(new_llh - old_llh) < 0.1)) break;
      old_llh = new_llh;
      // M-STEP: update the emission rates...
      MutSpecEM::get_all_map_spectra(with_ps);
    }
    if (loud ==1) cout<<endl;
    // test...
    if( gsl_matrix_isnonneg(spectra) == 0 || gsl_matrix_isnonneg(weights) == 0){
      printf("Error in MutSpecEM::run_em() after exact EM!\n");
      exit(1);
    }
  }  
  // Report results...
  double K = (double) Nsp*(Nch-1);
  double n = (double) Nsa;
  rllh = new_llh;
  rbic = (double) 2.0*new_llh - K*log(n);// BIC
  if (loud == 1){
    printf("BIC = %.6e\n", rbic);
    cout<<flush;
  }
  return(0);
}

// set the distance of the current reconstruction to the observation (X-x.mu^t)...
void MutSpecEM::set_dist(){
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, weights, spectra, 0.0, dist);
  gsl_matrix_mul_elements( dist, opp);
  gsl_matrix_add( dist, counts);
}

// Mean error in the reconstruction...
double MutSpecEM::get_mean_dev(){
  MutSpecEM::set_dist();
  double mean_dev = 0;
  for (int m=0; m<Nsa; m++){
    for(int j=0; j<Nch; j++){
      mean_dev += pow( gsl_matrix_get( dist, m, j), 2);
    }
  }
  mean_dev /= (Nch*Nsa);
  mean_dev = sqrt(mean_dev);
  return( mean_dev);
}

 // get the initial guess for the mixtures...
int MutSpecEM::get_weights_guess(){
  MutSpecEM::set_overlap();
  int m;
  int err=0;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic,2) default(shared)
#endif
  for (m=0; m<Nsa; m++){
    gsl_permutation * p = gsl_permutation_alloc(Nsp);
    int sign;
    gsl_matrix * L1  = gsl_matrix_alloc( Nsp, Nch);
    gsl_matrix * LHM = gsl_matrix_alloc( Nsp, Nsp);
    gsl_vector * cts  = gsl_vector_calloc(Nch);
    gsl_vector * weights_guess = gsl_vector_alloc(Nsp);
    gsl_vector * RHS  = gsl_vector_calloc(Nsp);
    gsl_vector_view col;
    gsl_vector_view oppV = gsl_matrix_row( opp, m);
    if ( gsl_blas_dasum(&oppV.vector)==0){// no opportunity - no activity!
      for (int a=0; a<Nsp; a++){
	gsl_matrix_set( weights, m, a, 0.0);
      }
    }
    // no observations - MAP value due to pseudo-counts only!
    else if ( gsl_vector_get(tot_counts,m) == 0){
      double val = 0;
      for (int a=0; a<Nsp; a++) val += gsl_matrix_get(overlap,m,a);
      for (int a=0; a<Nsp; a++) gsl_matrix_set( weights, m, a, (double) Nsp / val );
      //for (int a=0; a<Nsp; a++) gsl_matrix_set( weights, m, a, 0.0 );
    }
    else{ // In all other cases, do this...
      //#pragma omp critical
      {
	gsl_matrix_get_row( cts, counts, m);
	gsl_matrix_memcpy( L1, spectra);
      }
      for (int j=0; j<Nch; j++){
	col = gsl_matrix_column( L1, j);
	gsl_vector_scale( &col.vector, gsl_matrix_get(opp,m,j));
      }
      gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, L1, L1, 0.0, LHM);
      gsl_linalg_LU_decomp( LHM, p, &sign);
      if (gsl_linalg_LU_det (LHM, sign) == 0){
	//cout<<"ERROR in get_weights_guess(): LHM matrix is singular!\n";
	//exit(1);
	err=1;
      }
      else{
	gsl_vector_mul( cts, &oppV.vector);
	// x^m is the solution to:   LHM . x = RHS
	gsl_blas_dgemv( CblasNoTrans, 1.0, spectra, cts, 0.0, RHS);
	gsl_linalg_LU_solve( LHM, p, RHS, weights_guess);  
	for (int a=0; a<Nsp; a++){
	  if ( gsl_vector_get(weights_guess,a) < 0) gsl_vector_set(weights_guess,a,0);// weights_min;
	}
	gsl_matrix_set_row( weights, m, weights_guess);
      }
    }
    // cleanup...
    gsl_matrix_free(LHM);
    gsl_matrix_free(L1);
    gsl_permutation_free(p);
    gsl_vector_free(weights_guess);
    gsl_vector_free(RHS);
    gsl_vector_free(cts);
  }// END PARALLEL FOR
  //
  // add small constant to be safe...
  weights_min = get_pos_min(weights);
  gsl_matrix_add_constant( weights, 1.0e-6 * weights_min);
  return(err);
}


// update the emission spectra...
int MutSpecEM::get_spectra_guess(){
  int j;
  int err=0;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic,2) default(shared)
#endif
  for (j=0; j<Nch; j++){    
    gsl_permutation * p = gsl_permutation_alloc(Nsp);;
    int sign; 
    gsl_matrix * L1  = gsl_matrix_alloc(Nsa,Nsp);
    gsl_matrix * LHM = gsl_matrix_alloc(Nsp,Nsp);
    gsl_vector * cts = gsl_vector_calloc(Nsa);
    gsl_vector * RHS = gsl_vector_calloc(Nsp);
    gsl_vector * spec_guess = gsl_vector_alloc(Nsp);
    gsl_vector_view row;
    //#pragma omp critical
    {
      gsl_matrix_memcpy( L1, weights);
    }
    for (int m=0; m<Nsa; m++){
      row = gsl_matrix_row(L1,m);
      gsl_vector_scale( &row.vector, gsl_matrix_get(opp,m,j) );
    }
    gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, L1, L1, 0.0, LHM);
    gsl_linalg_LU_decomp( LHM, p, &sign);
    if (gsl_linalg_LU_det (LHM, sign) == 0){
      //cout<<"ERROR in get_spectra_guess(): LHM matrix is singular!\n";
      //exit(1);
      err=1;
      //for (int a=0; a<Nsp; a++) gsl_matrix_set(spectra,a,j,1.0/96.0);
    }
    else{
      gsl_matrix_get_col( cts, counts, j);
      gsl_blas_dgemv( CblasTrans, 1.0, L1, cts, 0.0, RHS);
      gsl_linalg_LU_solve( LHM, p, RHS, spec_guess);
      for (int a=0; a<Nsp; a++){
	if ( gsl_vector_get(spec_guess,a) < 0) gsl_vector_set(spec_guess, a, 0);//spectra_min;
      }
      gsl_matrix_set_col( spectra, j, spec_guess);   
    }
    //cleanup  
    gsl_vector_free(spec_guess);
    gsl_vector_free(RHS);
    gsl_vector_free(cts);
    gsl_matrix_free(LHM);
    gsl_matrix_free(L1);
    gsl_permutation_free(p); 
  }// END PARALLEL FOR
  //
  // add small constant...
  spectra_min = get_pos_min(spectra);
  gsl_matrix_add_constant( spectra, 1.0e-6*spectra_min);
    // normalize the spectra...
  gsl_vector_view spec;
  double val;
  for (int a=0; a<Nsp; a++){
    spec = gsl_matrix_row( spectra, a);
    val=gsl_blas_dasum( &spec.vector);
    if (val<=0){
      cout<<"ERROR in get_spectra_guess(): val = "<<val<<endl;
      exit(1);
    }
    gsl_vector_scale( &spec.vector, 1.0 / val);
  }  
  return(err);
}

// get ALL the maximum a posteriori (MAP) weights...
void MutSpecEM::get_all_map_weights(int loud, int with_ps){
  MutSpecEM::set_overlap();
  if (with_ps == 1){
    MutSpecEM::set_pseudo_counts();
  }  
  if ( gsl_matrix_isnonneg(weights) == 0){
    cout<<"ERROR in get_all_map_weights()!"<<endl;
    exit(1);
  }  
  // now loop over samples to find saddle points by numerical optimization...
  gsl_vector_set_zero(fexp);
  int m;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif  // START PARALLEL FOR
  for (m=0; m<Nsa; m++){
    if (loud == 1 && (0 == m % 100 || m==Nsa-1))
      printf("\rMAP in sample %i of %i...", m+1, Nsa);
    cout<<flush;
    gsl_vector_view oppV = gsl_matrix_row(opp,m);
    if (gsl_blas_dasum(&oppV.vector) == 0){
      for (int a=0; a<Nsp; a++) gsl_matrix_set( weights, m, a, 0.0);
      gsl_vector_set(done_with_ps,m,0);
    }
    else if ( gsl_vector_get(tot_counts,m) == 0){
      double val=0;
      if (with_ps==1){
	for (int a=0; a<Nsp; a++) val += gsl_matrix_get(overlap,m,a);
	for (int a=0; a<Nsp; a++) gsl_matrix_set( weights, m, a, 1.0/val);
	gsl_vector_set(done_with_ps, m, 1);
      }
      else{
	for (int a=0; a<Nsp; a++) gsl_matrix_set( weights, m, a, 0);
	gsl_vector_set(done_with_ps, m, 0);
      }
    }
    else{
      Fpar localFpar;
      gsl_matrix * M = gsl_matrix_alloc(Nsp,Nch);
      for (int a=0; a<Nsp; a++){
	for (int j=0; j<Nch; j++){
	  gsl_matrix_set( M, a, j, gsl_matrix_get(spectra,a,j)*gsl_matrix_get(opp,m,j));
	}
      }
      localFpar.M = M;
      localFpar.with_ps = with_ps;
      gsl_vector_view cts = gsl_matrix_row( counts, m);
      localFpar.cts = &cts.vector;
      gsl_vector * wts = gsl_vector_calloc(Nsp);
      gsl_matrix_get_row( wts, weights, m);
      gsl_vector * ps_cts = gsl_vector_calloc(Nsp);
      // with pseudo-counts ?!
      if ( with_ps==1 ){
	for (int a=0; a<Nsp; a++){
	  gsl_vector_set( ps_cts, a, gsl_matrix_get(pseudo_counts_E, m, a));
	}
	localFpar.ps_cts  = ps_cts;
	gsl_vector_set(done_with_ps, m, 1);
      }
      else{
	gsl_vector_set(done_with_ps, m, 0);
      }
      // get the MAP saddle point
      double val = get_map( wts, &localFpar);
      gsl_matrix_set_row( weights, m, wts);
      gsl_vector_set( fexp, m, val);
      gsl_vector_free(wts);
      gsl_vector_free(ps_cts);
      gsl_matrix_free(M);
    }
  }// END PARALLEL FOR
  //
  // test...
  if ( gsl_matrix_isnonneg(weights) == 0){
    cout<<"The new map weights have negative entries!\n";
    exit(1);
  }
}


// update the mutation spectra with their MAP estimate...
void MutSpecEM::get_all_map_spectra(int with_ps){
  int j;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic,2) default(shared)
#endif
  for (j=0; j<Nch; j++){// START PARALLEL FOR
    gsl_matrix * M = gsl_matrix_alloc(Nsp,Nsa);
    for (int a=0; a<Nsp; a++){
      for (int m=0; m<Nsa; m++){
	gsl_matrix_set( M, a, m, gsl_matrix_get(weights,m,a) * gsl_matrix_get(opp,m,j));
      }
    }
    Fpar localFpar;
    localFpar.M = M;
    localFpar.with_ps = with_ps;
    gsl_vector * ps_cts = gsl_vector_calloc(Nsp);
    if ( with_ps==1 ){
      for (int a=0; a<Nsp; a++){
	gsl_vector_set(ps_cts, a, gsl_matrix_get(pseudo_counts_M,j,a));
      }
      localFpar.ps_cts  = ps_cts;
    }
    // get all the counts in that channel...
    gsl_vector_view cts = gsl_matrix_column( counts, j);
    localFpar.cts = &cts.vector;
    // get the current spectra in that channel as start position...
    gsl_vector * spec = gsl_vector_alloc(Nsp);
    gsl_matrix_get_col( spec, spectra, j);
    get_map( spec, &localFpar);
    gsl_matrix_set_col( spectra, j, spec);
    gsl_matrix_free(M);
    gsl_vector_free(spec);
    gsl_vector_free(ps_cts);
  }// END PARALLEL FOR
  // test...
  if ( gsl_matrix_isnonneg(spectra) == 0){
    cout<<"New spectra are negative!"<<endl;
    exit(1);
  }
  // normalize the spectra...
  gsl_vector_view spectrum;
  for (int a=0; a<Nsp; a++){
    spectrum = gsl_matrix_row( spectra, a);
    gsl_vector_scale( &spectrum.vector, 1.0 / gsl_blas_dasum( &spectrum.vector));
  }  
}


// The actual numerical minimization routine to get the MAP estimates...
double get_map( gsl_vector * guess, Fpar * myFpar){
  // Here starts the minimizing step...
  int iter = 0, max_iter = 1.0e3; // max no. iterations
  int status;
  const gsl_multimin_fdfminimizer_type * T;
  gsl_multimin_fdfminimizer * s;
  gsl_multimin_function_fdf my_func;
  int nvar = (int) guess->size;
  // Give the elements of the function-to-minimize object...
  my_func.n      = nvar;
  my_func.f      = &F_map;
  my_func.df     = &dF_map;
  my_func.fdf    = &FdF_map;
  my_func.params = static_cast<void*>(myFpar);
  // Define type of minimization procedure...
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  //T = gsl_multimin_fdfminimizer_conjugate_pr; // to be used if above not available
  s = gsl_multimin_fdfminimizer_alloc( T, nvar);
  // Set and initialize the minimizer with LOG-TRANSFORMED DATA...
  gsl_vector * x  = gsl_vector_alloc(nvar);
  for (int a=0; a<nvar; a++){
    gsl_vector_set( x, a, log( gsl_vector_get(guess,a) ));
  }
  gsl_multimin_fdfminimizer_set( s, &my_func, x, 0.1, 0.01);
  // Now iterate to find the minimum...
  do{
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);     
    if (status) break;
    status = gsl_multimin_test_gradient( gsl_multimin_fdfminimizer_gradient(s), 1.0e-4);
  } 
  while (status == GSL_CONTINUE && iter < max_iter);
  
  double val = gsl_multimin_fdfminimizer_minimum(s);
  // copy the final result back into the argument, after EXPONENTIATION
  for (int a=0; a<nvar; a++) {
    gsl_vector_set( guess, a, exp( gsl_vector_get(gsl_multimin_fdfminimizer_x(s), a)));
  }
  // cleanup...
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
  return(val);
}



// the function to be minimized for finding the MAP estimates (x -> exp(x))
double F_map (const gsl_vector * x, void * p){
  Fpar * param = static_cast<Fpar*> (p);
  int N = (int) (param->cts)->size;
  int n = (int) x->size;
  double Fval = 0.0, arg=0.0;
  for (int j=0; j<N; j++){
    arg = 0.0;
    for  (int a=0; a<n; a++){
      arg += gsl_matrix_get(param->M, a, j) * exp( gsl_vector_get(x,a) );
    }
    if (arg>0) Fval += arg - gsl_vector_get(param->cts, j) * log(arg);
  }
  if (param->with_ps == 1){
    for  (int a=0; a < n; a++){
      Fval -= gsl_vector_get(param->ps_cts, a) * gsl_vector_get(x,a);
    }
  }
  return(Fval);
}



// this sets the gradient of above function (again in log-transformed coordinates)...
void dF_map (const gsl_vector * x, void * p, gsl_vector * grad){
  Fpar * param = static_cast<Fpar*> (p);
  int N = (int) (param->cts)->size;
  int n = (int) x->size;
  double val=0,sum;
  double * lambda = new double [N];
  double * expx = new double [n];
  for (int a=0;a<n;a++) expx[a] = exp( gsl_vector_get(x,a) );
  // a-independent term...
  for (int j=0; j<N; j++){
    lambda[j] = 0;
    for (int a=0; a<n; a++){
      lambda[j] += expx[a] * gsl_matrix_get( param->M, a, j);
    }
  }
  // now set the elements of the gradient...  
  gsl_vector_set_zero(grad);
  for (int a=0; a<n; a++){
    sum=0;
    for (int j=0; j<N; j++){
      val = gsl_matrix_get( param->M, a, j);
      if (  lambda[j] > 0 ) val *= (1.0 - gsl_vector_get( param->cts, j) / lambda[j] );
      sum += val;
    } 
    grad->data[a] = sum * expx[a];
  }
  if (param->with_ps == 1){
    for (int a=0; a<n; a++){
      gsl_vector_set( grad, a, gsl_vector_get(grad,a) - gsl_vector_get(param->ps_cts,a));
    }
  }
  delete [] lambda;
  delete [] expx;
}

// optimized variant of the two above...
void FdF_map (const gsl_vector * x, void * p, double * val, gsl_vector * grad){
  Fpar * param = static_cast<Fpar*> (p);
  int N = (int) (param->cts)->size;
  int n = (int) x->size;
  double * lambda = new double [N];
  double * expx = new double [n];
  double var=0,sum;
  for (int a=0;a<n;a++) expx[a] = exp( gsl_vector_get(x,a) );
  double Fval=0;
  // a-independent term...
  for (int j=0; j<N; j++){
    lambda[j] = 0;
    for (int a=0; a<n; a++){
      lambda[j] += expx[a] * gsl_matrix_get( param->M, a, j);
    }
    // now get the function value...
    if (lambda[j]>0) Fval += lambda[j] - gsl_vector_get(param->cts,j) * log(lambda[j]);
  }
  // now set the elements of the gradient...  
  gsl_vector_set_zero(grad);
  for (int a=0; a<n ; a++){
    sum = 0;
    for (int j=0; j<N; j++){
      var = gsl_matrix_get( param->M, a, j);
      if (lambda[j] > 0 ) var *= (1.0 -  gsl_vector_get(param->cts,j) / lambda[j]);
      sum += var;
    }
    grad->data[a] = sum * expx[a];
  }
  if (param->with_ps == 1){
    for (int a=0; a<n; a++){
      gsl_vector_set( grad, a, gsl_vector_get(grad,a) - gsl_vector_get(param->ps_cts,a));
      Fval -= gsl_vector_get(param->ps_cts,a) * gsl_vector_get(x,a);
    }
  }
  *val = Fval;
  delete [] lambda;
  delete [] expx;
}

	       

// the log-likelihood via saddle point approximation...
double MutSpecEM::get_llhood(){
  double llhood = 0, norm, val;
  gsl_permutation * p = gsl_permutation_alloc(Nsp);
  int sgn = 1;
  gsl_matrix * Hess = gsl_matrix_calloc( Nsp, Nsp);
  gsl_matrix * HessM = gsl_matrix_calloc( Nsp, Nsp);
  gsl_matrix * specM = gsl_matrix_calloc( Nsp, Nch);
  gsl_vector_view  oppV;//col
  // The llh is additive over independent samples...
  for (int m=0; m<Nsa; m++){
    oppV = gsl_matrix_row(opp,m);
    if (gsl_blas_dasum(&oppV.vector)==0) continue;
    // Saddle Point Approximation to the log-likelihood...
    llhood += (double) Nsp * log(2*PI) / 2;
    // Subtract fexp, i.e. the function values at the maximum..
    llhood -= gsl_vector_get( fexp, m);
    // constants...
    llhood -= gsl_vector_get( factorials, m);
    llhood += gsl_vector_get( ps_const, m);
    // now calculate the Hessian of f...
    gsl_matrix_set_zero(Hess);
    gsl_matrix_memcpy( specM, spectra);
    for (int j=0; j<Nch; j++){
      norm = 0;
      for (int a=0; a<Nsp; a++){
	norm += gsl_matrix_get( weights, m, a) * gsl_matrix_get( spectra, a, j);
      }
      if (norm == 0){
	printf("In get_llhood(): norm is %e!\n",norm);
	exit(1);
      }
      for (int a=0; a<Nsp; a++){
	for (int b=0; b<Nsp; b++){
	  val  = gsl_matrix_get(spectra,a,j);
	  val *= gsl_matrix_get(spectra,b,j);
	  val *= gsl_matrix_get(counts,m,j);
	  val /= pow(norm,2);
	  gsl_matrix_set( Hess, a, b, gsl_matrix_get(Hess,a,b) + val);
	}	
      }
    }
    // if this was arrived at with pseudo counts, then add this term... 
    if (gsl_vector_get(done_with_ps,m)==1){
      for (int a=0; a<Nsp; a++){
	gsl_matrix_set( Hess, a, a, 
			gsl_matrix_get(Hess,a,a) 
			+ gsl_matrix_get(pseudo_counts_E, m, a)
			/ pow( gsl_matrix_get(weights,m,a), 2)
			);
      }
    }
    if (gsl_matrix_isnonneg(Hess) == 0){
      exit(1);
    }
    gsl_matrix_memcpy( HessM, Hess);
    double Hmax = gsl_matrix_max(Hess);
    gsl_matrix_scale( Hess, 1.0 / Hmax);
    // Now compute the LU-decomposition and the log(|det(Hess)|)
    gsl_linalg_LU_decomp( Hess, p, &sgn);
    if( gsl_linalg_LU_det(Hess, sgn) == 0){
      printf("\nIn get_llhood for m=%i: det(Hf) = %e, ln|det(Hf)| = %e\nHf:\n", 
	     m+1, gsl_linalg_LU_det(Hess, sgn), gsl_linalg_LU_lndet(Hess));
      for (int a=0; a<Nsp; a++){
	for (int b=0; b<Nsp; b++){
	  printf("%.10e ", gsl_matrix_get(HessM,a,b) );
	}
	cout<<endl;
      }
      exit(1);
    }
    else{
      val = gsl_linalg_LU_lndet( Hess) + Nsp * log(Hmax);
      llhood -= 0.5*val ;    
    }
  }
  // cleanup...
  gsl_permutation_free(p);
  gsl_matrix_free(Hess);
  gsl_matrix_free(HessM);
  gsl_matrix_free(specM);
  // return
  return(llhood);
}

// use a local Gaussian approximation of the data-likelihood to obtain analytical error estimates...
void  MutSpecEM::get_error_estimates(int with_ps){
  if (with_ps == 1){
    MutSpecEM::set_pseudo_counts();
  } 
  gsl_matrix * Hessian = gsl_matrix_calloc(Nsp*Nch, Nsp*Nch);
  gsl_permutation * p = gsl_permutation_alloc(Nsp*Nch);
  int sgn = 1;
  double val, x, ol;
  gsl_vector_view wts, spec;
  // fill the entries of the (large) spectra Hessian matrix...
  for (int a=0; a<Nsp; a++){
    for (int b=0; b<Nsp; b++){
      for (int j=0; j<Nch; j++){
	val=0;
	spec = gsl_matrix_column(spectra,j);
	for (int m=0; m<Nsa; m++){
	  x = gsl_matrix_get(counts,m,j) * gsl_matrix_get(weights,m,a) * gsl_matrix_get(weights,m,b);
	  wts = gsl_matrix_row(weights,m);
	  gsl_blas_ddot( &spec.vector, &wts.vector, &ol);
	  x /= pow(ol,2); 
	  val += x;
	  if (with_ps && a==b){
	    x = gsl_matrix_get( pseudo_counts, m, a*Nch + j);
	    x /= pow(gsl_matrix_get(spectra,a,j),2);
	    val += x;
	  }
	}
	gsl_matrix_set(Hessian, a*Nch+j, b*Nch+j, val);
      }
    }
  }
  if (gsl_matrix_isnonneg(Hessian) == 0){
    cout<<"ERROR in get_error_estimates(): found negative value in the Hessian\n";
    exit(1);
  }
  // Now do LU-decomposition...
  gsl_linalg_LU_decomp( Hessian, p, &sgn);
  // and now fill the elements of the spectra error estimate...
  gsl_vector * inv_col  = gsl_vector_alloc(Nsp*Nch);
  gsl_vector * unit_vec = gsl_vector_alloc(Nsp*Nch);
  for (int a=0; a<Nsp; a++){
    for (int j=0; j<Nch; j++){
      gsl_vector_set_zero(unit_vec);
      gsl_vector_set(unit_vec, a*Nch+j, 1.0);
      gsl_linalg_LU_solve( Hessian, p, unit_vec, inv_col);  
      x = gsl_vector_get(inv_col,a*Nch+j);
      gsl_matrix_set( spectra_err, a,j, sqrt(x));
    }
  }
  // cleanup...
  gsl_permutation_free(p);
  gsl_matrix_free(Hessian);
  gsl_vector_free(unit_vec);
  gsl_vector_free(inv_col);
}





// Random step on a simplex, using a uniform jump distribution
void random_step_angle_uniform( gsl_vector * angle, double eps){
  int n = angle->size;
  double p, val;
  for( int i=0; i<n; i++){
    p = (double) rand() / RAND_MAX;
    p = (1.0 - 2.0*p) * eps;
    val = gsl_vector_get(angle,i) + p;
    if (val < 0){
      val = -val;
    }
    else if (val > PI/2){
      val = PI - val;
    }
    gsl_vector_set( angle, i, val);
  }
}


// Transform spherical angle coordinates to simplex coordinates...
double angle_to_simplex( gsl_vector * angle, gsl_vector * simplex){
  int n = angle->size;
  if (n+1 != (int) simplex->size){
    cout<<"ERROR in angle_to_simplex(): dimesions of two arguments incompatible!"<<endl;
    exit(1);
  }
  double * sinV = new double[n];
  double * cosV = new double[n];
  for ( int i =0; i<n; i++){
    sinV[i] = sin(gsl_vector_get(angle,i));
    cosV[i] = cos(gsl_vector_get(angle,i));
  }
  double mem=1;
  for ( int i =0; i<n; i++){
    gsl_vector_set( simplex, i, mem * cosV[i]);
    mem = mem*sinV[i];
  }
  gsl_vector_set( simplex, n, mem);
  for ( int i=0; i<=n; i++){
    mem = gsl_vector_get(simplex,i);
    gsl_vector_set( simplex, i, mem*mem);
  }
  double norm = gsl_blas_dasum(simplex);
  if (norm <=0){
    cout<<"norm = "<<norm<<endl;
    exit(1);
  }
  gsl_vector_scale( simplex, 1.0 / norm);
  double dE = 0.0;
  for ( int i=0; i<n; i++){
    dE += log(cosV[i]) + double(2*n - 2*i + 1) * log(sinV[i]);
  }
  dE += double(n+1)*log(2.0);
  delete [] sinV;
  delete [] cosV;
  return(dE);
}


// Transform simplex coordinates to angle coordinates ...
void simplex_to_angle( gsl_vector * simplex, gsl_vector * angle){
  int N = angle->size + 1;
  if (N != (int) simplex->size){
    cout<<"ERROR in simplex_to_angle(): dimesions of two arguments incompatible!"<<endl;
    exit(1);
  }
  double partSum = gsl_vector_get( simplex, N-1) + gsl_vector_get( simplex, N-2);
  double val=0;
  double x = gsl_vector_get( simplex, N-1);
  if (x > 0.0){
    val =  (sqrt(partSum) + sqrt(gsl_vector_get( simplex, N-2))) / sqrt(x);
    gsl_vector_set( angle, N-2,  PI - 2.0*atan(val));
  }
  else if (x == 0.0){
    gsl_vector_set( angle, N-2, 0.0);
  }
  else{
    cout<<"ERROR in simplex_to_angle()\n";
    exit(1);
  }
  //
  for ( int i=3; i<N+1; i++){
    val = sqrt( gsl_vector_get(simplex,N-i) / partSum);
    val = PI/2.0 - atan(val);
    if (val<0||val>PI/2.0){
      cout<<"ERROR in simplex_to_angle: angle is "<<val<<"!\n";
      exit(1);
    }
    gsl_vector_set( angle, N-i, val);
    partSum += gsl_vector_get( simplex, N-i);
  }
}

// make a random step on the simplex...
void make_step_on_simplex( gsl_vector * simplex, double eps){
  int N = simplex->size;
  gsl_vector * angle = gsl_vector_alloc(N-1);
  gsl_vector * trial = gsl_vector_alloc(N-1);
  simplex_to_angle( simplex, angle);
  double oldE = angle_to_simplex(angle,simplex);
  double newE = 0;
  int acc=0;
  while (acc==0){
    gsl_vector_memcpy( trial, angle);
    random_step_angle_uniform( trial, eps);
    newE = angle_to_simplex(trial,simplex);
    if ( newE > oldE || (double) (rand()/RAND_MAX) < oldE-newE){
      acc=1;
    }
  }  
  gsl_vector_free(angle);
  gsl_vector_free(trial);
}


// MCMC in the space of mutational spectra (Nsp Nch-dimensional simplices)...
void MutSpecEM::run_mcmc(int mcmc_steps, const char * pre){
  // int mcmc_steps=1.0e5;
  double acc_rate=0, eps=1.0e-6, corr=0, scale = 1.0, delta=0;
  int acc=0, rej=0, samples=0;
  char buff[256];
  sprintf( buff, "%s_%i_MCMC_LLH_traj.txt", pre, Nsp);
  FILE * Etraj_fp = fopen( buff,"w");
  sprintf( buff,"%s_%i_MCMC_spectra_traj.txt", pre, Nsp);
  FILE * mutraj_fp = fopen(buff,"w");
  gsl_matrix * angles_trial    = gsl_matrix_alloc(Nsp,Nch-1);
  gsl_matrix * angles_current  = gsl_matrix_alloc(Nsp,Nch-1);
  gsl_matrix * spectra_current = gsl_matrix_alloc(Nsp,Nch);
  gsl_matrix * spectra_trial   = gsl_matrix_alloc(Nsp,Nch);
  gsl_matrix_set_zero(spectra_mean);
  gsl_matrix_set_zero(spectra_err);
  gsl_matrix * spectra_square = gsl_matrix_calloc(Nsp,Nch);
  gsl_matrix * mem = gsl_matrix_calloc(Nsp,Nch);
  gsl_vector * angle  = gsl_vector_alloc(Nch-1);
  gsl_vector * spectrum  = gsl_vector_alloc(Nch);
  for (int a=0; a<Nsp; a++){
    gsl_matrix_get_row( spectrum, spectra, a);
    simplex_to_angle( spectrum, angle);
    gsl_matrix_set_row( angles_current, a, angle);
  }
  gsl_matrix_memcpy( spectra_current, spectra);
  int with_ps=1;
  corr=0;
  for(int a=0; a<Nsp;a++){
    gsl_matrix_get_row( angle, angles_current, a);
    corr += angle_to_simplex(angle,spectrum);
  }
  MutSpecEM::get_all_map_weights( 0, with_ps);
  double oldE = MutSpecEM::get_llhood() + corr;
  double newE=0;//, dev=0;
  double p,crit;
  srand(time(NULL));
  printf("Start MCMC at E = %.6e\n", oldE);
  for (int step=0; step<mcmc_steps; step++){
    corr=0;
    // step in spectra space...
    for(int a=0; a<Nsp;a++){
      gsl_matrix_get_row( angle, angles_current, a);
      random_step_angle_uniform( angle, eps);
      gsl_matrix_set_row( angles_trial, a, angle);
      corr += angle_to_simplex(angle,spectrum);
      gsl_matrix_set_row( spectra_trial, a, spectrum);
    }
    // MutSpecEM::get_all_map_weights( 0, with_ps);
    // compute the data log-likelihood...
    //newE = MutSpecEM::get_llhood() + corr;
    newE = MutSpecEM::get_app_llhood( spectra_trial ) + corr;
    if (newE > oldE){
      crit = 1.0;
      p = 0.0;
    }
    else{
      crit = exp((newE - oldE));
      p = (double) rand() / RAND_MAX;
    }
    // test for acceptance...
    if ( p  < crit ){
      gsl_matrix_memcpy( angles_current, angles_trial);
      gsl_matrix_memcpy( spectra_current, spectra_trial);
      oldE = newE;
      acc++;
    }
    else{
      rej++;
    }
    acc_rate = (double) acc / (acc+rej);
    // update the reference spectrum and saddle point...
    if(acc % 10 == 0){
      gsl_matrix_memcpy( spectra, spectra_current);
      MutSpecEM::get_all_map_weights( 0, with_ps);
    }
    // report...
    if ((step+1) % 1000 == 0){
      printf("\rMCMC step %6i of %6i: E = %.6e, eps = %.3e, acc.rate = %.3f",
	     step+1,mcmc_steps, oldE, eps, acc_rate);
      cout<<flush;
    }
    // reset eps to achieve an accpetance rate of ca. 0.25...
    if (step % 1000 == 0){
      if (acc_rate < 0.25){
	if (delta>0){
	  scale *= 0.9;
	  delta=0;
	}
	eps = eps/(1.0+scale);
	delta--;
      }
      else{
	if (delta<0){
	  scale *= 0.9;
	  delta=0;
	}
	eps = eps*(1.0+scale);
	if (eps>1.0e-1) eps = 1.0e-1;
	delta++;
      }
      if (scale<0.2) scale = 0.2;
      acc=0;
      rej=0;
    }
    // record mean and second moment...
    if (step % 100 == 0){
      gsl_matrix_add( spectra_mean, spectra_current);
      gsl_matrix_memcpy( mem, spectra_current);
      gsl_matrix_mul_elements(mem,mem);
      gsl_matrix_add( spectra_square, mem);
      samples++;
      // print trajectory in spectra space...
      fprintf(Etraj_fp,"%i %.3e\n",step, oldE);
      for (int a=0;a<Nsp;a++){
	for (int j=0; j<Nch; j++){
	  fprintf(mutraj_fp,"%.3e ", gsl_matrix_get( spectra_current, a, j));
	}
      }
      fprintf(mutraj_fp,"\n");
    }
  }
  gsl_matrix_get_row( spectrum, spectra_mean, 0);
  double norm =  gsl_blas_dasum(spectrum);
  gsl_matrix_scale(spectra_mean, 1.0 / norm);
  gsl_matrix_memcpy(spectra,spectra_mean);
  MutSpecEM::get_all_map_weights( 0, with_ps);
  gsl_matrix_scale(spectra_square, 1.0 / norm);
  gsl_matrix_memcpy(mem,spectra_mean);
  gsl_matrix_mul_elements(mem,mem);
  gsl_matrix_scale(mem,-1.0);
  gsl_matrix_add(spectra_square,mem);
  for (int a=0; a<Nsp; a++){
    for (int j=0; j<Nch; j++){
      gsl_matrix_set(spectra_err,a,j, sqrt(gsl_matrix_get(spectra_square,a,j)));
    }
  }
  cout<<endl;
  //cleanup...
  fclose(Etraj_fp);
  fclose(mutraj_fp);
  gsl_matrix_free(angles_trial);
  gsl_matrix_free(angles_current);
  gsl_matrix_free(spectra_current);
  gsl_matrix_free(spectra_trial);
  gsl_matrix_free(spectra_square);
  gsl_matrix_free(mem);
  gsl_vector_free(angle);
  gsl_vector_free(spectrum);
}



// MCMC in the space of mutational spectra (Nsp Nch-dimensional simplices)...
void MutSpecEM::freeze_out(int steps){
  double acc_rate=0, eps=1.0e-6, corr=0, scale = 1.0, delta=0;
  int acc=0, rej=0;
  gsl_matrix * angles_trial    = gsl_matrix_alloc(Nsp,Nch-1);
  gsl_matrix * angles_current  = gsl_matrix_alloc(Nsp,Nch-1);
  gsl_matrix * spectra_current = gsl_matrix_calloc(Nsp,Nch);
  gsl_matrix * spectra_trial   = gsl_matrix_calloc(Nsp,Nch);
  gsl_vector * angle  = gsl_vector_alloc(Nch-1);
  gsl_vector * spectrum  = gsl_vector_alloc(Nch);
  for (int a=0; a<Nsp; a++){
    gsl_matrix_get_row( spectrum, spectra, a);
    simplex_to_angle( spectrum, angle);
    gsl_matrix_set_row( angles_current, a, angle);
  }
  gsl_matrix_memcpy( spectra_current, spectra);
  int with_ps=1;
  corr=0;
  for(int a=0; a<Nsp;a++){
    gsl_matrix_get_row( angle, angles_current, a);
    corr += angle_to_simplex(angle,spectrum);
  }
  MutSpecEM::get_all_map_weights( 0, with_ps);
  double oldE = MutSpecEM::get_llhood();// + corr;
  double newE=0;
  srand(time(NULL));
  printf("Start freeze out at E = %.6e\n", oldE);
  for (int step=0; step<steps; step++){
    corr=0;
    for(int a=0; a<Nsp;a++){
      gsl_matrix_get_row( angle, angles_current, a);
      random_step_angle_uniform( angle, eps);
      gsl_matrix_set_row( angles_trial, a, angle);
      corr += angle_to_simplex(angle,spectrum);
      gsl_matrix_set_row( spectra_trial, a, spectrum);
    }
    //MutSpecEM::get_all_map_weights( 0, with_ps);
    // compute the data log-likelihood...
    newE = MutSpecEM::get_app_llhood( spectra_trial );// + corr;
    //newE = MutSpecEM::get_llhood();// + corr;
    // test for acceptance...
    if (newE > oldE){
      gsl_matrix_memcpy( angles_current, angles_trial);
      gsl_matrix_memcpy( spectra_current, spectra_trial);
      // reset reference...
      gsl_matrix_memcpy( spectra, spectra_current);
      MutSpecEM::get_all_map_weights( 0, with_ps);
      oldE = newE;
      acc++;
    }
    else{
      rej++;
    }
    // acceptance rate...
    acc_rate = (double) acc / (acc+rej);
    // report...
    if ((step+1) % 1000 == 0){
      printf("\rfreeze out step %6i of %6i: E = %.6e, eps = %.3e, acc.rate = %.3f",
	     step+1, steps, oldE, eps, acc_rate);
      cout<<flush;
    }
    // reset eps to achieve an accpetance rate of ca. 0.25...
    if (step % 1000 == 0){
      if (acc_rate < 0.25){
	if (delta>0){
	  scale *= 0.9;
	  delta=0;
	}
	eps = eps/(1.0+scale);
	delta--;
      }
      else{
	if (delta<0){
	  scale *= 0.9;
	  delta=0;
	}
	eps = eps*(1.0+scale);
	if (eps>1.0e-1) eps = 1.0e-1;
	delta++;
      }
      if (scale<0.2) scale = 0.2;
      acc=0;
      rej=0;
    }
    if (eps < 1.0e-10) break;
  }
  gsl_matrix_memcpy(spectra,spectra_current);
  MutSpecEM::get_all_map_weights( 0, with_ps);
  cout<<endl;
  //cleanup...
  gsl_matrix_free(angles_trial);
  gsl_matrix_free(angles_current);
  gsl_matrix_free(spectra_current);
  gsl_matrix_free(spectra_trial);
  gsl_vector_free(angle);
  gsl_vector_free(spectrum);
}




// get the min positive element of a matrix...
double get_pos_min(gsl_matrix * M){
  double posmin=0,val;
  for (int i=0; i<(int) M->size1; i++){
    for (int j=0; j<(int) M->size2; j++){
      val = gsl_matrix_get(M,i,j);
      if (val > 0 && (val<posmin || posmin == 0)) posmin = val;
    }
  }
  return(posmin);
}
// ... or vector...
double get_pos_min(gsl_vector * v){
  double posmin=0,val;
  for (int i=0; i<(int) v->size; i++){
    val = gsl_vector_get(v,i);
    if (val > 0 && (val<posmin || posmin == 0)) posmin = val;
  }
  return(posmin);
}



void print_vector(const gsl_vector * x){
  for (int i=0; i<(int)x->size; i++){
    printf("%.2e ", gsl_vector_get(x,i));
  }
  cout<<endl;  
}
void print_matrix(const gsl_matrix * M){
  for (int i=0; i<(int)M->size1; i++){
    for (int j=0; j<(int)M->size2; j++){
      printf("%.2e ", gsl_matrix_get(M,i,j));
    }
    cout<<endl;
  }
}


// the log-likelihood via saddle point approximation...
double MutSpecEM::get_app_llhood(gsl_matrix * spectra_trial){
  double llhood = 0;
  if(gsl_matrix_ispos(weights) == 0){
    cout<<"ERROR in get_app_llhood(): weights is not positive definite matrix.\n";
    exit(1);
  }
  if(gsl_matrix_ispos(spectra_trial) == 0){
    cout<<"ERROR in get_app_llhood(): spectra is not positive definite matrix.\n";
    exit(1);
  }
  int m=0;
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic, 1) default(shared)
#endif
  for (m=0; m<Nsa; m++){//START PARALLEL FOR
    double lh=0;
    // allocate local variables...
    double norm, val;
    gsl_permutation * p = gsl_permutation_alloc(Nsp);
    int sgn = 1;
    gsl_matrix * pseudo = gsl_matrix_calloc( Nsp, Nch);
    gsl_vector * pseudo_sum = gsl_vector_calloc( Nsp);
    gsl_matrix * Hess = gsl_matrix_calloc( Nsp, Nsp);
    gsl_matrix * HessM = gsl_matrix_calloc( Nsp, Nsp);
    gsl_vector * grad = gsl_vector_alloc(Nsp);
    gsl_vector * wts_spec = gsl_vector_alloc(Nch);
    gsl_vector * mem = gsl_vector_alloc(Nsp);
    gsl_vector * y   = gsl_vector_alloc(Nsp);
    gsl_vector_view  oppV, wtsV, colV, ctsV;
    // start calculation...
    oppV = gsl_matrix_row(opp,m);
    wtsV = gsl_matrix_row(weights,m);
    if (gsl_blas_dasum(&oppV.vector)==0) continue;
    for (int j=0; j<Nch; j++){ 
      colV = gsl_matrix_column( spectra_trial, j);
      gsl_blas_ddot( &wtsV.vector, &colV.vector, &val);
      if (val>0){
	gsl_vector_set( wts_spec, j, val);
      }
      else{
	cout<<"ERROR in get_app_llhood(): wts_spec("<<j<<") = "<<val<<endl;
	exit(1);
      }
    }
    // Saddle Point Approximation to the log-likelihood...
    lh += (double) Nsp * log(2.0*PI) / 2.0;
    // Subtract the function value at the saddle point..
    gsl_blas_dgemv( CblasNoTrans, 1.0, spectra_trial, &oppV.vector, 0.0, mem);
    gsl_vector_memcpy( grad, mem);
    gsl_blas_ddot( &wtsV.vector, mem, &val);
    lh -= val;
    // pseudo-counts term...
    MutSpecEM::get_naive_ps_cts( pseudo, spectra_trial, &oppV.vector);
    for (int a=0; a<Nsp; a++){
      ctsV = gsl_matrix_row(pseudo,a);
      norm = gsl_blas_dasum(&ctsV.vector);
      gsl_vector_set( pseudo_sum, a, norm);
    } 
    double incr = 0;
    for (int j=0; j<Nch; j++){ 
      if (gsl_matrix_get(counts,m,j)>0){
	incr +=  gsl_matrix_get(counts,m,j) * log( gsl_vector_get(wts_spec,j) * gsl_matrix_get(opp,m,j) );
      }
      for (int a=0; a<Nsp; a++){
	val  = gsl_matrix_get( weights, m, a);
	val *= gsl_matrix_get( spectra_trial, a, j);
	val *= gsl_matrix_get( opp, m, j);
	if (gsl_matrix_get(pseudo,a,j) > 0 && val > 0){
	  incr += gsl_matrix_get(pseudo,a,j) * log(val); 
	}
      }
    }
    lh += incr;
    // factorials...
    lh -= gsl_vector_get( factorials, m);
    // now calculate the gradient of f...
    for (int a=0; a<Nsp; a++){
      val = gsl_vector_get(grad,a);
      for (int j=0; j<Nch; j++){ 
	val -= gsl_matrix_get(counts,m,j) * gsl_matrix_get(spectra_trial,a,j) / gsl_vector_get(wts_spec,j);
      }
      val -= gsl_vector_get(pseudo_sum,a) / gsl_matrix_get(weights,m,a);
      gsl_vector_set(grad,a,val);
    }
    // now calculate the Hessian matrix of f...
    gsl_matrix_set_zero(Hess);
    for (int j=0; j<Nch; j++){
      for (int a=0; a<Nsp; a++){
	for (int b=0; b<Nsp; b++){
	  gsl_matrix_set(HessM, a, b, gsl_matrix_get(spectra_trial,a,j) * gsl_matrix_get(spectra_trial,b,j));
	}
      }
      val = gsl_matrix_get(counts,m,j);
      val /= pow( gsl_vector_get(wts_spec,j), 2);
      gsl_matrix_scale( HessM, val);
      gsl_matrix_add( Hess, HessM);
    }    
    // pseudo-counts...
    for (int a=0; a<Nsp; a++){
      val  = gsl_matrix_get(Hess,a,a);
      val += gsl_vector_get(pseudo_sum,a) / pow( gsl_matrix_get(weights,m,a), 2); 
      gsl_matrix_set(Hess,a,a,val);
    }
    //gsl_matrix_memcpy( HessM, Hess);
    double Hmax = gsl_matrix_max(Hess);
    if (Hmax <= 0){
      cout<<"ERROR in get_app_llhood(): Hmax = "<<Hmax<<endl;
      exit(1);
    }
    gsl_matrix_scale( Hess, 1.0 / Hmax);
    // Now compute the LU-decomposition and the log(|det(Hess)|)
    gsl_linalg_LU_decomp( Hess, p, &sgn);
    if( gsl_linalg_LU_det(Hess, sgn) == 0){
      cout<<"ERROR in get_app_llhood(): det(Hess) = 0\n";
      exit(1);
    }
    else{
      // Hessian term: log det(Hess)
      val = (double) (gsl_linalg_LU_lndet( Hess) + Nsp * log(Hmax));
      lh -= 0.5*val;
      // gradient term... 
      // using: log(grad^T Hess^{-1} grad) = log(y^T Hess^T y), with Hess y = grad
      gsl_linalg_LU_solve( Hess, p, grad, y);
      gsl_vector_memcpy(grad,y);
      gsl_blas_dgemv(CblasTrans,1.0,Hess,y,0.0,mem);
      gsl_blas_ddot( grad, mem, &val);
      lh += 0.5*val / Hmax;
    }
    // cleanup local variables...
    gsl_permutation_free(p);
    gsl_matrix_free(Hess);
    gsl_matrix_free(pseudo);
    gsl_matrix_free(HessM);
    gsl_vector_free(pseudo_sum);
    gsl_vector_free(wts_spec);
    gsl_vector_free(grad);
    gsl_vector_free(mem);
    gsl_vector_free(y);
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      llhood += lh;
    }
  }  
  // return
  return(llhood);
}



void MutSpecEM::get_naive_ps_cts( gsl_matrix * Pseudo, gsl_matrix * Spectra, gsl_vector * Opp){
  double norm = 0;
  for (int a=0; a<Nsp; a++){
    for (int j=0; j<Nch; j++){
      gsl_matrix_set( Pseudo, a, j, gsl_matrix_get(Spectra,a,j) * gsl_vector_get(Opp,j));
      norm += gsl_matrix_get(Pseudo,a,j);
    }
  }
  if (norm>0){
    gsl_matrix_scale(Pseudo, (double) Nsp / norm);
  }
  else{
    cout<<"ERROR in get_naive_ps_cts(): norm = "<<norm<<endl;
    exit(1);
  }
}
