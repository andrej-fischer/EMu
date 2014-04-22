/* 

******************************************************************************

Copyright (c) 09-16-13  Genome Research Ltd.

Author: Andrej Fischer (af7[at]sanger.ac.uk)

This file is part of EMu

EMu is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation; either version 3 
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  
If not, see <http://www.gnu.org/licenses/>.

*/


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <ctype.h> 
#include <string>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif


// GSL headers...
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

using namespace std;

struct cmdl_opts{
  const char * mut_file_name; 
  const char * cnv_file_name; 
  const char * reg_file_name; 
  const char * pre;
  const char * suff;
  const char * chr_dir;
  string bin_str;
  int bin_size, max_chr, no_chr, cnv_def;
};

// get command line arguments...
void get_opts( int argc, const char ** argv, cmdl_opts& opts);

// mapping mutations to 96-channels
void mut_to_chn( cmdl_opts& opts,
		 map<string,char>& mut_idx, 
		 map<char,char>& base_idx);

// mutational opportunity per fixed-sized bin

void get_opp( map<char,char>& base_idx, 
	      cmdl_opts& opts);

// collect mutations per fixed-size bin
void get_mut(  cmdl_opts& opts,
	       map<string,char>& mut_idx, 
	       map<char,char>& base_idx
	       );

// get file handles to all the chromosomes
void get_chrs(ifstream* ifs, int * start, int * stride, cmdl_opts& opts);

// *** MAIN START ***
int main (int argc, const char * argv[]){
  cmdl_opts opts;
  // Read options from the command line...
  get_opts( argc, argv, opts);
  map<char,char> base_idx;
  base_idx['N'] =  0;
  base_idx['A'] =  1;
  base_idx['C'] =  2;
  base_idx['G'] =  3;
  base_idx['T'] =  4;
  // mutation channel numbering...
  map<string, char> mut_idx;
  string key [] = {"C>A","G>T","C>G","G>C","C>T","G>A","T>A","A>T","T>C","A>G","T>G","A>C"};
  char   val [] = {0,0,1,1,2,2,3,3,4,4,5,5};
  for (int i=0; i<12; i++) mut_idx.insert(pair<string,char>(key[i],val[i]));
  // Get the opportunity per sample and genomic segment including copy number information...
  if (opts.chr_dir != NULL && opts.cnv_file_name != NULL){
    get_opp( base_idx, opts);
  }
  // Now fetch the locations and channels of observed mutations
  // check whether the mutations are already in a 96-channel format...
  if (opts.mut_file_name != NULL){
    cout<<"Translating mutations to 96 chn format..."<<flush;
    get_mut( opts, mut_idx, base_idx);
    cout<<"done"<<endl;
  }
  return(0);
}
//*** MAIN END ***


void get_opts( int argc, const char ** argv, cmdl_opts& opts){
  opts.reg_file_name = NULL;
  opts.mut_file_name = NULL;
  opts.cnv_file_name = NULL;
  opts.chr_dir       = NULL;
  opts.pre           = "out";
  opts.suff          = NULL;
  opts.bin_size      = 0;
  opts.max_chr       = 100;
  opts.no_chr        = 0;
  opts.cnv_def       = 2; //copy number default
  int opt_idx = 1;
  string opt_switch;
  while ( opt_idx < argc && (argv[opt_idx][0] == '-')){
    opt_switch = argv[opt_idx];
    if ( opt_switch.compare("--regions") == 0){
      opt_idx++;
      opts.reg_file_name = argv[opt_idx];
    }
    else if ( opt_switch.compare("--mut") == 0){
      opt_idx++;
      opts.mut_file_name = argv[opt_idx];
    }
    else if ( opt_switch.compare("--chr") == 0){
      opt_idx++;
      opts.chr_dir = argv[opt_idx];
    }
    else if ( opt_switch.compare("--cnv") == 0){
      opt_idx++;
      opts.cnv_file_name = argv[opt_idx];
    } 
    else if ( opt_switch.compare("--pre") == 0){
      opt_idx++;
      opts.pre = argv[opt_idx];
    } 
    else if ( opt_switch.compare("--suff") == 0){
      opt_idx++;
      opts.suff = argv[opt_idx];
    } 
    else if ( opt_switch.compare("--bin") == 0){
      opt_idx++;
      opts.bin_size = int(atof( argv[opt_idx]));
    }
    else if ( opt_switch.compare("--default") == 0){
      opt_idx++;
      opts.cnv_def = atoi( argv[opt_idx]);
    }
    else {
      cout<<"ERROR: unexpected option "<< opt_switch <<endl;
      cout << "For usage, see README" << endl;
      exit(1);
    }
    opt_switch.clear();
    opt_idx++;
  }
  // *** TESTS ***
  if (opts.pre == NULL && opts.bin_size>0){
    cout<<"ERROR: you have to set --pre [prefix] for the bin-wise output files!\n";
    exit(1);
  }
  if (opts.cnv_file_name == NULL && opts.chr_dir == NULL){
    cout<<"ERROR: for mutational opportunity, you have to set --cnv [cnv_file] and --chr [dir]\n";
    exit(1);
  }
  if (opts.reg_file_name != NULL && opts.bin_size > 0 ){
    cout<<"ERROR: you have to decide whether to get data per bin genome wide (with --bin [size]) ";
    cout<<"OR whether to get data collapsed over specific regions (with --regions [reg_file])\n";
    exit(1);
  }
  if (opts.bin_size>0){
    char buff[128];
    sprintf( buff, "per-%ib", opts.bin_size);
    opts.bin_str.assign(buff);
    int bs = opts.bin_size;
    if (bs == 1000)     opts.bin_str = "per-1kb";
    if (bs == 10000)    opts.bin_str = "per-10kb";
    if (bs == 100000)   opts.bin_str = "per-100kb";
    if (bs == 1000000)  opts.bin_str = "per-1Mb";
    if (bs == 10000000) opts.bin_str = "per-10Mb";
  }
}


void get_chrs(ifstream* ifs, int * start, int * stride, cmdl_opts& opts){
  if (opts.chr_dir == NULL){
    cout<<"ERROR: you must set the directory of chromosome sequences with --chr [dir].\n";
    exit(1);
  }
  char chr_fn [128];
  //char c;
  string line;
  if (opts.no_chr > 0){
    for (int i=0; i<opts.no_chr; i++){
      ifs[i].close();
    }
    opts.no_chr = 0;
  }  
  for (int i=1; i<=opts.max_chr; i++){
    sprintf( chr_fn, "%s/chr%i.fa", opts.chr_dir, i);  
    ifs[i-1].open(chr_fn,ios::in);
    if ( ifs[i-1].fail() ){
      ifs[i-1].close();
      break;
    }
    opts.no_chr++;
    getline( ifs[i-1], line);
    if (line[0] != '>'){
      printf("ERROR: first line in %s should start with '>'.\n", chr_fn);
      exit(1);
    }
    start[i-1]  = (int) ifs[i-1].tellg();
    stride[i-1] = 0;
    getline( ifs[i-1], line);
    stride[i-1] = (int) line.length();
    ifs[i-1].seekg(0,ifs[i-1].beg);
  }
}



void get_regions( vector<int> * region_start,
		  vector<int> * region_stop,
		  cmdl_opts& opts
		  ){
  ifstream reg_ifs;
  reg_ifs.open(opts.reg_file_name,ios::in);
  int reg_start,reg_stop;
  string line;
  stringstream line_ss;
  char chr_buff[8];
  while( (line.clear(), getline( reg_ifs, line)) ){
    line_ss.clear();
    line_ss.str(line);
    // get chromosome
    line_ss.getline( chr_buff, 8, ' ');
    int chr = atoi(chr_buff) - 1;
    // get postion
    line_ss >> reg_start >> reg_stop;
    // check
    if (region_start[chr].size() > 0 && reg_start < region_start[chr].back()){
      printf("ERROR: regions in %s need to be ordered by start coordinate in chromosome\n", opts.reg_file_name);
      cout<<line<<endl;
      printf("%i %i\n", region_start[chr].back(),region_stop[chr].back());
      exit(1);
    }
    if (region_start[chr].size() > 0 && reg_start <= region_stop[chr].back()){
      reg_start = region_stop[chr].back() + 1;
    }
    if (reg_start > reg_stop){
      continue;
    }
    region_start[chr].push_back(reg_start);
    region_stop[chr].push_back(reg_stop);
  }
  reg_ifs.close();
}


//*******************************************************************************************
// translate mutation information to the 96 tri-nucleotide channel format...
void get_mut( cmdl_opts& opts,
	      map<string,char>& mut_idx, 
	      map<char,char>& base_idx
	      ){
  int no_chn = 96;
  ifstream * chr_ifs = new ifstream [opts.max_chr];
  int * start        = new int [opts.max_chr];
  int * stride       = new int [opts.max_chr];
  get_chrs(chr_ifs,start,stride,opts);
  int coord=0, idx=0, chr=0, pos=0;
  char left, middle, right, contxt=0;
  ifstream mut_ifs;
  mut_ifs.open( opts.mut_file_name, ios::in);
  string mut_chn_fn(opts.mut_file_name);
  mut_chn_fn.append(".96");
  FILE * mut_chn_fp = fopen( mut_chn_fn.c_str(), "w" );
  string line, mut, sample;
  istringstream line_ss;
  char wt;
  map<string,int*> mut_per_sample; 
  map<string,gsl_matrix**> mut_per_bin;
  // *** GET REGION INFORMATION ***
  vector<int> * region_start = NULL;
  vector<int> * region_stop  = NULL;
  if (opts.reg_file_name != NULL){
    region_start = new vector<int> [opts.no_chr];
    region_stop  = new vector<int> [opts.no_chr];
    get_regions(region_start,region_stop,opts);
  }
  int last=0,no_char=0,no_bins=0,no_bases=0,chn=0;
  unsigned int curr_bin=0;
  char nucl [] = {'a','c','g','t'};
  char Nucl [] = {'A','C','G','T','a','c','g','t'};
  std::set<char> nt(nucl,nucl+4);
  std::set<char> NT(Nucl,Nucl+8);
  // *** MUTATIONS ***
  //expects format 
  // sample chr coord A>G
  while( (line.clear(), getline( mut_ifs, line)) ){
    if (line.empty()) break;
    // read the mutation information, expects the format: "1234 A>G"
    line_ss.clear();
    line_ss.str(line);
    line_ss >> sample;
    // *** FIRST ENCOUNTER OF THIS SAMPLE ***
    if (mut_per_sample.count(sample) == 0){
      int * muts = new int [no_chn];
      for (int j=0; j<no_chn; j++) muts[j] = 0;
      mut_per_sample.insert(pair<string,int*>(sample,muts));
      if (opts.bin_size>0){
	gsl_matrix ** mut_pb = new gsl_matrix * [opts.no_chr];
	for (int chr=0; chr<opts.no_chr; chr++){
	  // compute the number of bins in this chromosome...
	  chr_ifs[chr].seekg( 0, chr_ifs[chr].end);
	  last     = chr_ifs[chr].tellg();
	  no_char  = last - start[chr] + 1;
	  no_bases = no_char - int((double) no_char / (stride[chr]+1));	
	  no_bins  = (int) ceil( (double) no_bases / opts.bin_size);
	  mut_pb[chr] = gsl_matrix_calloc( no_bins, no_chn);
	}
	mut_per_bin.insert(pair<string,gsl_matrix**>(sample,mut_pb));
      }
    }
    // get chromosome
    line_ss >> chr;
    chr--;
    if (start[chr] == 0 || stride[chr] == 0){
      printf("ERROR: no information available for chr %i\n%s\n", chr+1,line.c_str());
      exit(1);
    }
    // get chromosome coordinate...
    line_ss >> coord;
    // *** test if coord falls into one of the regions ***
    if (opts.reg_file_name != NULL){
      int found=0;
      for (int reg_ct=0; reg_ct < (int) region_start[chr].size(); reg_ct++){
	if (coord >= region_start[chr][reg_ct] && coord <= region_stop[chr][reg_ct]){
	  found = 1;
	  break;
	}
      }
      if (found == 0) continue;
    }
    // get the actual mutation
    line_ss >> mut;
    if (coord <= 0) {
      printf("\nERROR: Expect 1-based coordinate system: %i %s (?=%s)\n",
	     coord, mut.c_str(), line.c_str());
      exit(1);
    }
    // get the position of the base BEFORE the one in question in the file...
    if (coord<2) continue;
    // coord-2 is the 0-based coordinate of the base 5' of the mutated base
    pos = start[chr] + int( (double) (coord-2) / stride[chr])*(stride[chr]+1) + ((coord-2) % stride[chr]);
    (chr_ifs[chr]).seekg( pos, (chr_ifs[chr]).beg);
    (chr_ifs[chr]).get(left);
    if (left == '\n') (chr_ifs[chr]).get(left);
    (chr_ifs[chr]).get(middle);
    if (middle == '\n') (chr_ifs[chr]).get(middle);
    (chr_ifs[chr]).get(right);
    if (right == '\n') (chr_ifs[chr]).get(right);
    if (nt.count(left))   left   -= 32; // convert to upper case
    if (nt.count(middle)) middle -= 32; 
    if (nt.count(right))  right  -= 32;
    //
    //if (middle != 'A' && middle != 'C' && middle != 'G' && middle != 'T'){
    if ( NT.count(middle) == 0 ){
      printf("\nERROR: Mutation %s at position %i in chr %i in %s not possible! Reference is '%c'.\n",
	     mut.c_str(), coord, chr+1, sample.c_str(), middle);
      line.clear();
      continue;
    }
    if (left == 'N' || right == 'N') continue;
    //
    wt = mut[0];// the wild type base (for 'A>T' this would be 'A')
    if ( base_idx[wt] != base_idx[middle] && base_idx[wt] != (5-base_idx[middle]) ){
      printf("ERROR: Initial base (%c) at %i in chr %i different from reference (%c) in %s!\n", 
	     mut[0], coord, chr+1, middle, sample.c_str());
      line.clear();
      continue;
    }
    // deciding the context...
    if (middle == 'C' || middle == 'T'){
      contxt = (base_idx[left]-1)*4 + (base_idx[right]-1);
    }
    else if (middle == 'G' || middle == 'A'){
      contxt = (5 - base_idx[right] - 1)*4 + (5 - base_idx[left] - 1);
    }
    else{
      cout<<"ERROR:1\n";
      exit(1);
    }
    // get the actual mutation channel...
    idx = mut_idx[mut];
    chn = idx*16 + contxt;
    fprintf( mut_chn_fp, "%s %i %i %s %i\n", 
	     sample.c_str(), chr+1, coord, mut.c_str(), chn);
    mut_per_sample[sample][chn]++;
    if (opts.bin_size>0){
      curr_bin = (unsigned int) floor((double) (coord-1) / opts.bin_size);
      double val = gsl_matrix_get( mut_per_bin[sample][chr], curr_bin, chn);
      gsl_matrix_set( mut_per_bin[sample][chr], curr_bin, chn, val + 1);
    }
    line.clear();
  }
  // cleanup...
  fclose( mut_chn_fp);
  mut_ifs.close();
  // print collected mutations...
  string samples_fn(opts.mut_file_name);
  samples_fn.append(".samples");
  FILE * samples_fp = fopen( samples_fn.c_str(), "w" );
  string mat_fn(opts.mut_file_name);
  mat_fn.append(".mut.matrix");
  FILE * mat_fp = fopen( mat_fn.c_str(), "w" );
  map<string,int*>::iterator it;
  for (it = mut_per_sample.begin(); it != mut_per_sample.end(); ++it){
    fprintf(samples_fp, "%s\n", (it->first).c_str());
    for (int j=0; j<no_chn; j++){
      fprintf(mat_fp, "%i ", (it->second)[j]);
    }
    fprintf(mat_fp, "\n");
  }
  fclose(mat_fp);
  fclose(samples_fp);
  // *** print mutations per bin ***
  if (opts.bin_size > 0){
    map<string,gsl_matrix**>::iterator it;
    for (it = mut_per_bin.begin(); it != mut_per_bin.end(); ++it){
      for (int chr=0; chr<opts.no_chr; chr++){
	string out_fn(opts.pre);
	out_fn.append(".");
	out_fn.append(it->first);
	out_fn.append(".mut.");
	char chrbuff [32];
	sprintf( chrbuff, "chr%i.%s.txt", chr+1, opts.bin_str.c_str());
	out_fn.append(chrbuff);
	FILE * out_fp = fopen(out_fn.c_str(),"w");
	for (int bin=0; bin<(int) (it->second)[chr]->size1; bin++){
	  for (int j=0; j<no_chn; j++){
	    fprintf( out_fp, "%i ", (int) gsl_matrix_get( (it->second)[chr], bin, j));
	  }
	  fprintf( out_fp, "\n");
	}
	fclose(out_fp);
	gsl_matrix_free((it->second)[chr]);
      }
    }
  }
  else if (opts.reg_file_name != NULL){
    for (int chr=0; chr<opts.no_chr; chr++){
      region_start[chr].clear();
      region_stop[chr].clear();
    }
    delete [] region_start;
    delete [] region_stop;
  }
}


//*******************************************************************************************
// translate sequences to opportunity tracks, i.e. which mutational channels are open...
void get_opp( map<char,char>& base_idx, 
	      cmdl_opts& opts
	      ){
  //
  int no_chn = 96;
  // get chrs file handles
  ifstream * chr_ifs = new ifstream [opts.max_chr];
  int * start        = new int [opts.max_chr];
  int * stride       = new int [opts.max_chr];
  for (int i=0; i<opts.max_chr; i++){
    start[i]  = 0;
    stride[i] = 0;
  }
  get_chrs( chr_ifs, start, stride, opts);
  // copy number registers
  map<string,vector<int>*> cnv_start;
  map<string,vector<int>*> cnv_stop;
  map<string,vector<int>*> cnv_mult;
  vector<string> samples;
  ifstream cnv_ifs;
  cnv_ifs.open( opts.cnv_file_name, ios::in);    
  string line,sample;
  stringstream line_ss;
  char buff[128];
  int cstart,cstop,cmult,chr;
  // get cnv information
  map<string,string> sample_fn;
  int no_samples=0;
  // *** GET COPY NUMBER INFORMATION ***
  while( (line.clear(), getline( cnv_ifs, line)) ){
    if (line.length() == 0) break;
    line_ss.clear();
    line_ss.str(line);
    // get sample name
    line_ss >> buff;//white space separated
    sample.assign(buff);
    if ( cnv_start.count(sample) == 0){
      no_samples++;
      samples.push_back(sample);
      vector<int> * starts = new vector<int> [opts.no_chr];
      vector<int> * stops  = new vector<int> [opts.no_chr];
      vector<int> * mults  = new vector<int> [opts.no_chr];
      cnv_start.insert(pair<string,vector<int>*>(sample,starts));
      cnv_stop.insert(pair<string,vector<int>*>(sample,stops));
      cnv_mult.insert(pair<string,vector<int>*>(sample,mults));
      if (opts.pre != NULL){
	string ofn(opts.pre);
	ofn.append(".");
	ofn.append(sample);
	sample_fn.insert(pair<string,string>(sample,ofn));
      }
    }
    line_ss >> buff;
    chr = atoi(buff) - 1;
    if (start[chr] == 0){
      printf("ERROR: no information for chr %i: %s\n", chr+1, line.c_str());
      exit(1);
    }
    line_ss >> cstart >> cstop >> cmult;//one-based
    if (cnv_start[sample][chr].size() > 0){
      if (cstart <= cnv_stop[sample][chr].back() || cstop <= cstart){
	printf("ERROR: %s", line.c_str());
	exit(1);
      }
      if (cnv_stop[sample][chr].back() == -2){
	printf("Didn't expect another cnv region in sample %s chr %i\n", sample.c_str(), chr+1);
	exit(1);
      }
    }
    cnv_start[sample][chr].push_back(cstart-1);//zero-based
    cnv_stop[sample][chr].push_back(cstop-1);
    cnv_mult[sample][chr].push_back(cmult);//mulitplier
  }
  cnv_ifs.close();
  //
  // fall back to default if no explicit cnv regions given
  for (int s=0; s<no_samples;s++){
    string sample = samples[s];
    for(int chr=0; chr<opts.no_chr; chr++){
      if (cnv_start[sample][chr].empty()){
	cnv_start[sample][chr].push_back(0);
	cnv_stop[sample][chr].push_back(-2);
	cnv_mult[sample][chr].push_back(opts.cnv_def);//default multiplier
      }
      //
    }
  }
  //exit(0);
  // collapsed opportunity per sample: allocation
  unsigned int ** opp_per_sample = new unsigned int * [no_samples];
  for (int s=0; s<no_samples; s++){
    opp_per_sample[s] = new unsigned int [no_chn];
    for (int j=0; j<no_chn; j++) opp_per_sample[s][j] = 0;
  }
  // *** GET REGION INFORMATION ***
  vector<int> * region_start = NULL;
  vector<int> * region_stop  = NULL;
  if (opts.reg_file_name != NULL){
    region_start = new vector<int> [opts.no_chr];
    region_stop  = new vector<int> [opts.no_chr];
    get_regions(region_start,region_stop,opts);
  }
  //
  // *** CHROMOSOMES ***
#ifdef _OPENMP
#pragma omp parallel for schedule( dynamic,1) default(shared)
#endif
  for(chr = 0; chr < opts.no_chr; chr++){
    // skip if no regions in this chromosome
    if (opts.reg_file_name != NULL){
      if (region_start[chr].empty()) continue;
    }
    // opportunity in this chromosome...
    unsigned int ** opp_per_chr = new unsigned int * [no_samples];
    for (int s=0; s<no_samples; s++){
      opp_per_chr[s] = new unsigned int [no_chn];
      for (int j=0; j<no_chn; j++) opp_per_chr[s][j] = 0;
    }
    int last=0,no_char=0,no_bins=0,no_bases=0;
    unsigned int curr_bin=0;
    double val;
    int idx;
    char left, middle, right, contxt=0;
    // compute the number of bins in this chromosome...
    chr_ifs[chr].seekg( 0, chr_ifs[chr].end);
    last     = chr_ifs[chr].tellg();
    no_char  = last - start[chr];
    no_bases = no_char - int((double) no_char / (stride[chr]+1));//new-line chars don't count!
    //exit(0);
    gsl_matrix * OppTrack = NULL; 
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      printf("Computing all opportunity for chr %2i (%i bases)", chr+1, no_bases);
      if (opts.bin_size>0){
	no_bins  = (int) ceil( (double) no_bases / opts.bin_size);
	printf(" in %6i bins", no_bins);
	// allocate the opportunity track...
	OppTrack = gsl_matrix_calloc( no_samples, no_bins*no_chn);
      }
      cout<<endl; 
    } 
    int * next_start  = new int [no_samples];
    int * next_stop   = new int [no_samples];
    int * copy_number = new int [no_samples];
    list<int> loci_of_interest;
    for (int s=0; s<no_samples; s++){
      next_start[s] = 0;
      next_stop[s]  = 0;
      copy_number[s] = opts.cnv_def;//default copy number
      if ( (cnv_stop[samples[s]])[chr].back() == -2){
	(cnv_stop[samples[s]])[chr].pop_back();
	(cnv_stop[samples[s]])[chr].push_back(no_bases);
      }
      //dummy, never used
      (cnv_start[samples[s]])[chr].push_back(no_bases);
      (cnv_stop[samples[s]])[chr].push_back(no_bases);
      (cnv_mult[samples[s]])[chr].push_back(opts.cnv_def);
      // get loci of interest...
      for (int i=0; i< (int) (cnv_start[samples[s]])[chr].size(); i++){
	loci_of_interest.push_back((cnv_start[samples[s]])[chr][i]);
	loci_of_interest.push_back((cnv_stop[samples[s]])[chr][i] + 1);
      }
    } 
    loci_of_interest.sort();
    loci_of_interest.unique();
    int j0=0, old_bin=-1;
    left   = 'N';
    middle = 'N';
    right  = 'N';
    // go to start of the chromosome sequence...
    chr_ifs[chr].seekg( start[chr], chr_ifs[chr].beg);
    int next_reg = 0;
    int first = (opts.reg_file_name == NULL) ? 0 : max( 0, region_start[chr][0] - 3);
    last      = (opts.reg_file_name == NULL) ? no_bases : region_stop[chr].back() + 3;
    // *** LOCI ***
    for (int ct = first; ct < last; ct++){
      // This changes the copy number if needed...
      if (ct == loci_of_interest.front()){
	for (int s=0; s<no_samples; s++){
	  if( ct == (cnv_start[samples[s]])[chr][next_start[s]] ){
	    copy_number[s] = (cnv_mult[samples[s]])[chr][next_start[s]];
	    next_start[s]++;
	  }
	  else if (ct == (cnv_stop[samples[s]])[chr][next_stop[s]] + 1){
	    copy_number[s] = opts.cnv_def;
	    next_stop[s]++;
	  }
	}
	loci_of_interest.pop_front();
      }
      // test whether within one of the regions
      if (opts.reg_file_name != NULL){
	if (ct < region_start[chr][next_reg]){
	  continue;
	}
	else if (ct == region_start[chr][next_reg]){
	  left   = 'N';// reset at the beginning of a new region
	  middle = 'N';
	  right  = 'N';
	}
	else if (ct > region_stop[chr][next_reg]){
	  next_reg++;
	  continue;
	}
      }
      // get the next base...
      chr_ifs[chr].get(right); 
      if (right == '\n') chr_ifs[chr].get(right); 
      if (right!='A' && right!='C' && right!='G' && right!='T') right = 'N';
      // the tri-nucleotide must not contain a 'N'!
      if ( ct>0 && left != 'N' && middle != 'N' && right != 'N'){
	// deciding the context...
	if (middle == 'C' || middle == 'T'){
	  contxt = (base_idx[left]-1)*4 + (base_idx[right] - 1);
	}
	else if (middle == 'G' || middle == 'A'){// the pyrimidine (C,T) is on the reverse strand
	  contxt = (5 - base_idx[right] - 1)*4 + (5 - base_idx[left] - 1);// the contxt is also reversed
	}
	// deciding and updating the three channels...
	j0 = (middle == 'C' || middle == 'G') ? 0 : 3;
	if (OppTrack != NULL){
	  // e.g. for bins_size = 1000, ct = 0..999 will be in curr_bin = 0
	  curr_bin = (unsigned int) floor((double) (ct-1) / opts.bin_size);
	  if ((int) curr_bin > old_bin){
	    //printf("\rIn bin %6i...", curr_bin+1);
	    //cout<<flush;
	    old_bin=curr_bin;
	  }
	}
	for (int s=0; s < no_samples; s++){
	  for (int j = j0; j < j0+3; j++){
	    opp_per_chr[s][j*16+contxt] += copy_number[s];
	  }
	}
	if (OppTrack != NULL){
	  for ( int s=0; s < no_samples; s++){
	    for (int j = j0; j < j0+3; j++){
	      idx = curr_bin * no_chn + j*16 + contxt;
	      val = gsl_matrix_get( OppTrack, s, idx);
	      gsl_matrix_set( OppTrack, s, idx, val + copy_number[s]);
	    }
	  }
	}
      }
      left = middle;
      middle = right;
    } 
    // Now collect the opportunity in this chromosome...
#ifdef _OPENMP
#pragma omp critical
#endif
    {
      for (int s=0; s<no_samples; s++){
	for (int j=0; j<no_chn; j++){
	  opp_per_sample[s][j] += opp_per_chr[s][j];
	}
	delete [] opp_per_chr[s];
      }
      delete [] opp_per_chr;
    }
    // bin-wise
    if (OppTrack != NULL){
      for ( int s=0; s< no_samples; s++){
	string fn = sample_fn[samples[s]];
	char chrbuff [32];
	sprintf( chrbuff, ".opp.chr%i.%s.txt", chr+1, opts.bin_str.c_str());
	fn.append(chrbuff);
	FILE * track_fp = fopen( fn.c_str(), "w");
	// BINS:
	for (int i=0; i<no_bins; i++){
	  // CHANNELS:
	  for (int j=0; j<no_chn; j++){
	    fprintf( track_fp, "%i ", (int) gsl_matrix_get( OppTrack, s, i*no_chn + j));
	  }
	  fprintf( track_fp, "\n");
	}
	fclose(track_fp);
      }   
      gsl_matrix_free(OppTrack);
    }
    //local clean up  
    delete [] copy_number;
    delete [] next_start;
    delete [] next_stop;
    loci_of_interest.clear();
    //printf("done\n");
  }
  // print total opportunity in all chromosomes...
  string samples_fn(opts.cnv_file_name);
  samples_fn.append(".samples");
  FILE * samples_fp = fopen( samples_fn.c_str(), "w" );
  string mat_fn(opts.cnv_file_name);
  mat_fn.append(".opp.matrix");
  FILE * mat_fp = fopen( mat_fn.c_str(), "w" );
  for (int s=0; s<no_samples; s++){
    fprintf(samples_fp, "%s\n", samples[s].c_str());
    for (int j=0; j<no_chn; j++){
      fprintf(mat_fp, "%i ", opp_per_sample[s][j]);
    }
    fprintf(mat_fp, "\n");
    delete [] opp_per_sample[s];
  }
  delete [] opp_per_sample;
  fclose(mat_fp);
  fclose(samples_fp);
  cnv_start.clear();
  cnv_stop.clear();
  cnv_mult.clear();
}

