/* Copyright (C) 2020, Murad Banaji
 *
 * This file is part of COVIDAGENT v0.2
 *
 * COVIDAGENT is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 3, 
 * or (at your option) any later version.
 *
 * COVIDAGENT is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with COVIDAGENT: see the file COPYING.  If not, see 
 * <https://www.gnu.org/licenses/>

 */

 /* 
 * Some description of the modelling carried out
 * using COVIDAGENT can be found at 
 * maths.mdx.ac.uk/research/modelling-the-covid-19-pandemic/
 */


#include "inf.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> // random seeding
#include <ctype.h>
#include <string.h>
#include <iostream>
#include <random>

// Maximum number of infected individuals - memory limit?
#define MAXINFS 10000000

#define max(A, B) ((A) > (B) ? (A) : (B))

//Externally declared (bad practice I know!)

int inflist[MAXINFS]; //For book-keeping free spaces in list
std::default_random_engine generator;

int getline(FILE *fp, char s[], int lim)
{
  /* store a line as a string, including the terminal newline character */
  int c=0, i;

  for(i=0; i<lim-1 && (c=getc(fp))!=EOF && c!='\n';++i)
    s[i] = c;
  if (c == '\n') {
    s[i] = c;
    ++i;
  }
  s[i] = '\0';
  return i;
}

//From https://stackoverflow.com/questions/29787310/does-pow-work-for-int-data-type-in-c
int int_pow(int base, int exp){
  int result = 1;
  while (exp)
    {
      if (exp % 2)
	result *= base;
      exp /= 2;
      base *= base;
    }
  return result;
}

int getnthblock(char *s, char *v, int len, int n){
  // get the nth valid block (only spaces count as separators) from a string s and put it in v. Returns the next position in  s. 
  int i, j, k;
  i=0, k=0;
  if(len < 2){
    fprintf(stderr, "ERROR in getnthblock in inf2.cc: third argument must be at least 2.\n");
    return 0;
  }
  v[0] = '\0'; // in case we have an empty string, return an empty word.
  while(s[k] != '\0'){
    j=0;
    while(isspace((int) s[k])) // skip space
      k++;
    for(i=0;i<n-1;i++){
      while(!(isspace((int) s[k]))) // skip first word 
        k++;
      while(isspace((int) s[k])) //skip space
        k++;
    }
    while((j<len-1) && !(isspace((int) s[k]))){ // get the word
      v[j++] = s[k++];
    }
    v[j++] = '\0';
    if(j==len){
      fprintf(stderr, "WARNING: in routine getnthblock in file inf2.cc: word is longer than maximum length.\n");
    }
    return k;
  }
  return 0;
}

FILE *openftowrite(const char fname[]){
  FILE *fd;
  if(!(fd=fopen(fname, "w"))){
    fprintf(stderr, "FILE \"%s\" could not be opened for writing. EXITING.\n", fname);exit(0);
  }
  return fd;
} 
FILE *openftoread(const char fname[]){
  FILE *fd;
  if(!(fd=fopen(fname, "r"))){
    fprintf(stderr, "FILE \"%s\" could not be opened for reading. EXITING.\n", fname);exit(0);
  }
  return fd;
} 


int readDataFile(const char fname[], int **data, int max){
  //No error checking
  FILE *fd;
  int lim=1000;
  char oneline[lim];
  char val1[50], val2[50], val3[50];
  int numlines=0, len;

  fd=openftoread(fname);
  while((len = getline(fd, oneline, lim)) > 0 && numlines<max){
    if ((oneline[0] == '#') || (oneline[0] == '/') || (oneline[0] == '\n') || (oneline[0] == '\0')){} // comment/empty lines
    else{
      getnthblock(oneline, val1, 50, 1);
      data[numlines][0]=atoi(val1);
      getnthblock(oneline, val2, 50, 2);
      data[numlines][1]=atoi(val2);
      getnthblock(oneline, val3, 50, 3);
      data[numlines][2]=atoi(val3);
      numlines++;
    }
  }
  if(numlines>=max){
    fprintf(stderr, "Data file \"%s\" too long to be read. EXITING.\n", fname);
    exit(0);
  }

  fclose(fd);

  return numlines;

}


int getoption(char *fname, const char optname[], int num, char v[], int max){
  FILE *fd;
  int len, j, flag=0;
  char oneline[200];
  char modname[50];
  int lim=200;
  fd = openftoread(fname);
  while((len = getline(fd, oneline, lim)) > 0){
    j=0;
    while((isspace((int) oneline[j])) || (oneline[j] == 13)){j++;}
    if ((oneline[j] == '/') || (oneline[j] == '\n') || (oneline[j] == '\0')){} // comment/empty lines
    else{
      getnthblock(oneline, modname, 50, 1);
      if(strcmp(modname, optname) == 0){
        flag = 1;
        getnthblock(oneline, modname, 50, num+1);
        if((int)(strlen(modname)) < max-1)
	  strcpy(v, modname);
        else{
	  fprintf(stderr, "ERROR in routine getoption: Option %s in file %s has value %s which is too long.\n", optname, fname, modname);
	  v[0] = '\0';
	  fclose(fd);
	  return -1;
        }
        break;
      }
    }
  }
  fclose(fd);
  if(flag==0){
    fprintf(stderr, "WARNING in routine getoption: Option %s could not be found in file %s. Setting to default value.\n", optname, fname);
    v[0] = '\0';
    return -2;
  }
  return 0;

}

int getoptioni(char *fname, const char optname[], int defval, FILE *fd1){
  char tempword[200];
  int val;
  if(getoption(fname, optname, 1, tempword, 200)!=0)
    val=defval;
  else
    val=atoi(tempword);
  fprintf(fd1, "#%s %d\n", optname, val);
  return val;
}

int getoption2i(char *fname, const char optname[], int defval, FILE *fd1){
  char tempword[200];
  int val;
  if(getoption(fname, optname, 2, tempword, 200)!=0)
    val=defval;
  else
    val=atoi(tempword);
  fprintf(fd1, "#%s %d\n", optname, val);
  return val;
}

float getoptionf(char *fname, const char optname[], float defval, FILE *fd1){
  char tempword[200];
  float val;
  if(getoption(fname, optname, 1, tempword, 200)!=0)
    val=defval;
  else
    val=atof(tempword);
  fprintf(fd1, "#%s %.4f\n", optname, val);
  return val;
}

float getoption2f(char *fname, const char optname[], float defval, FILE *fd1){
  char tempword[200];
  float val;
  if(getoption(fname, optname, 2, tempword, 200)!=0)
    val=defval;
  else
    val=atof(tempword);
  fprintf(fd1, "#%s %.4f\n", optname, val);
  return val;
}


int randnum(int max){
  return rand()%max;
}

int randpercentage(double perc){// to 1 d.p. Casting to int is flooring
  int intperc=(int)(10.0*perc);
  //fprintf(stderr, "%d\n", intperc);
  if (randnum(1000)<intperc)
    return 1;
  return 0;
}



long factorial(int x){
  int i;
  long factx = 1;
  for(i=1; i<=x ; i++ )
    factx *= i;
  return factx;
}

//Poisson distribution (see Wikipedia page)
double Pois(double lamb, int k){
  return exp(k*log(lamb) - lamb -lgamma(k+1));
}

//Geometric distribution including 0
double Geom(double Ep, int k){
  return pow((1.0-1.0/(1.0+Ep)), (double)k)/(1.0+Ep);
}

//binomial distribution shifted to centre at 0
double binom(int n, int param){
  if(param==0){//no distribution
    if(n==0)
      return 1.0;
    else
      return 0.0;
  }
  else if(param==2){//one each side
    if(n==-1 || n==1)
      return 0.25;
    else if(n==0)
      return 0.5;
    else
      return 0.0;
  }
  else if(param==4){//two each side
    if(n==-2 || n==2)
      return 0.0625;
    else if(n==-1 || n==1)
      return 0.25;
    else if(n==0)
      return 0.375;
    else
      return 0.0;
  }
  else if(param==6){//three each side
    if(n==-3 || n==3)
      return 0.015625;
    else if(n==-2 || n==2)
      return 0.09375;
    else if(n==-1 || n==1)
      return 0.234375;
    else if(n==0)
      return 0.3125;
    else
      return 0.0;
  }
  else{
    fprintf(stderr, "ERROR - invalid parameter in binom. EXITING.\n");
    exit(0);
  }
  return 0.0;
}

//Choose from binomial distribution (even parameter up to 6)
int choosefrombin(int param){
  int r;
  int tot;
  if(param==0)
    return 0;
  r=randnum(1000)+1; //1 to 1000
  tot=(int)(1000.0*binom(3,param));
  if(r<tot)
    return -3;
  tot+=(int)(1000.0*binom(3,param));
  if(r<tot)
    return 3;
  tot+=(int)(1000.0*binom(2,param));
  if(r<tot)
    return -2;
  tot+=(int)(1000.0*binom(2,param));
  if(r<tot)
    return 2;
  tot+=(int)(1000.0*binom(1,param));
  if(r<tot)
    return -1;
  tot+=(int)(1000.0*binom(1,param));
  if(r<tot)
    return 1;
  return 0;

}



//Poisson distribution cast to integers (not rounded) and with the tail given to largest value
int *discPois(double lamb, int max){
  int *P;
  int k;
  P=(int *)malloc((size_t) ((max+1)*sizeof(int)));
  P[0]=1000;
  for(k=1;k<=max;k++){
    P[k]=P[k-1]-(int)(1000.0*Pois(lamb, k-1));
  }
  return P;
}



//Geometric distribution cast to integers (not rounded) and with the tail given to zero
int *discGeom(double lamb, int max){
  int *P;
  int k;
  P=(int *)malloc((size_t) ((max+1)*sizeof(int)));
  P[0]=1000;
  for(k=1;k<=max;k++){
    P[k]=P[k-1]-(int)(1000.0*Geom(lamb, k-1));
  }
  return P;
}

//Sample from a distribution with (cast to integers out of 1000)
int choosefromdist(int P[], int totP){//P has totp+1 entries
  int r=randnum(1000)+1; //1 to 1000
  int ct=totP;
  while(ct>0){
    if(r<P[ct])
      return ct;
    ct--;
  }
  return 0;

}

//
// Gamma distribution
//

double gamma(double shp, double scl, std::default_random_engine & generator)
{
  static std::gamma_distribution<double> dist;
  return dist(generator, std::gamma_distribution<double>::param_type(shp, scl));
}

//
// Uniform distribution
//

double unif(double lend, double rend, std::default_random_engine & generator)
{
  static std::uniform_real_distribution<double> dist1;
  return dist1(generator, std::uniform_real_distribution<double>::param_type(lend, rend));
}

//
// uniform distribution on integers
//

int unifi(int lend, int rend, std::default_random_engine & generator)
{
  static std::uniform_int_distribution<int> dist1;
  return dist1(generator, std::uniform_int_distribution<int>::param_type(lend, rend));
}

//
// normal distribution
//

double norml(double mean, double stdev, std::default_random_engine & generator)
{
  static std::normal_distribution<double> dist1;
  return dist1(generator, std::normal_distribution<double>::param_type(mean, stdev));
}



// infected class
// inf P is the distribution, maxP is the greatest num to infect
// time steps in days
// infected at age 0, but updates to 1 at the start of each timestep
// quarantined? Then can't infect
// ill? Can die if ill
// How many will this person infect? Follows the distribution
// The assumption is that it isn't just randomly about contact
// perhaps also biology?
// The infection times are chosen from a uniform distribution on some
// range of values

//In case of any user-defined distribution
inf::inf(int orgnum, int P[], int maxP){
  age = 0;
  ill = 0;
  quar = 0;
  num = orgnum;
  numtoinf=choosefromdist(P, maxP);
}

//in case of gamma distribution
inf::inf(int orgnum, double shp, double scl){
  age = 0;
  ill = 0;
  quar = 0;
  num = orgnum;
  double number = gamma(shp, scl, generator);
  if((numtoinf=int(round(number)))>MAXDISCPROB-1)
    numtoinf=MAXDISCPROB-1;
//fprintf(stderr, "numtoinf[%d]=%d\n", orgnum, numtoinf);
}


// Set the times at which infection occurs: uniform distribution C++ generator
void inf::setinftimes(int rmin, int rmax){
  int i; 
  for(i=0;i<numtoinf;i++){
    inftimes[i]=unifi(rmin, rmax, generator);
    (infnums[inftimes[i]])++;
  }
}

//Set the times at which infection occurs: gamma distribution
void inf::setinftimes(double shp, double scl){
  int i; 
  int num;
  for(i=0;i<numtoinf;i++){
    num = int(round(gamma(shp, scl, generator)));
    if(num>0 && num<MAXAGE)
      inftimes[i]=num;
    else if(num>=MAXAGE)
      inftimes[i]=MAXAGE-1;
    else
      inftimes[i]=1;
    (infnums[inftimes[i]])++;
  }
}

// The first available position in the list
int firstfreepos(int *inflist, int maxinfs){
  int i;
  for(i=0;i<MAXINFS;i++){
    if(inflist[i]==0)
      return i;
  }
  return -1;
}

inf **infar(long nl, long nh)
/* allocates a set of pointers to infs */
{
  long nrow = nh-nl+1;
  inf **m;
  m=(inf **) malloc((size_t)((nrow+1)*sizeof(inf*)));
  if (!m) fprintf(stderr, "allocation failure in infar()\n");
  m+=1;
  m-=nl;

  return m;
}

void free_infar(inf **m, long nl, long nh)
/* free an inf array allocated by infar() */
{
  free((char *) (m+nl-1));
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+1)*sizeof(int*)));
  if (!m) fprintf(stderr, "allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+1)*sizeof(int)));
  if (!m[nrl]) fprintf(stderr, "allocation failure 2 in matrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch){
  free((char *) (m[nrl]+ncl-1));
  free((char *) (m+nrl-1));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+1)*sizeof(double*)));
  if (!m) fprintf(stderr, "allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+1)*sizeof(double)));
  if (!m[nrl]) fprintf(stderr, "allocation failure 2 in matrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch){
  free((char *) (m[nrl]+ncl-1));
  free((char *) (m+nrl-1));
}


int nextpos=0;// only needed if not freeing

int create(inf *infs[], int gamswtch, double alpha, double beta, int P[], int maxP, int inf_gam, int inf_start, int inf_end, double inf_mid, double inf_tm_shp, int *numinf, int *numcurinf, int *newinfs, int *numill, double percill, double percdeath, double time_to_death, double dist_on_death, double time_to_recovery, double dist_on_recovery, double time_to_sero, double dist_on_sero, double quardate, double quarp, double dist_on_quardate, double testp, double testdelay, double testdelay_shp){
  int i=firstfreepos(inflist, MAXINFS);
  //int i=nextpos++;
  int j;
  double inf_scl=inf_mid/inf_tm_shp;
  if(i==-1){//no more space
    fprintf(stderr, "Ran out of space in list - you can consider resetting MAXINFS. EXITING.\n");
    exit(0);
  }
  if(!gamswtch)
    infs[i] = new inf(i, P, maxP);
  else
    infs[i] = new inf(i, alpha, beta);
  (*numinf)++;(*numcurinf)++;(*newinfs)++;

  if(dist_on_sero>=0)//discrete simple
    infs[i]->sero_time=(int)time_to_sero+choosefrombin((int)dist_on_sero);
  else//normal dist., -dist_on_sero=stdev
    infs[i]->sero_time=int(round(norml(time_to_sero, -dist_on_sero, generator)));

  for(j=0;j<MAXAGE;j++){//number to infect at time j
    infs[i]->infnums[j]=0;
  }
  // who falls ill?
  // Currently unused - left in for potential use
  if(randpercentage(percill)){
    if(randpercentage(percdeath)){
      infs[i]->ill=-1;//falls ill and dies
      if(dist_on_death>=0)//discrete simple
	infs[i]->dth_time=(int)time_to_death+choosefrombin((int)dist_on_death);
      else//normally distributed, -dist_on_death=stdev
	infs[i]->dth_time=int(round(norml(time_to_death, -dist_on_death, generator)));

    }
    else{
      infs[i]->ill=1;//falls ill but recovers
      if(dist_on_recovery>=0)
	infs[i]->recov_time=(int)time_to_recovery+choosefrombin((int)dist_on_recovery);
      else{//normal dist, -dist_on_recovery=stdev
	infs[i]->recov_time=int(round(norml(time_to_recovery, -dist_on_recovery, generator)));
	// if(infs[i]->recov_time>MAXAGE)
	//   fprintf(stderr, "recov_time=%d\n", infs[i]->recov_time);
      }
    }
    (*numill)++;
    //fprintf(stderr, "ill=%d\n", infs[i]->ill);
  }
  else{//won't fall ill
    if(dist_on_recovery>=0)
      infs[i]->recov_time=(int)time_to_recovery+choosefrombin((int)dist_on_recovery);
    else//normal dist, -dist_on_recovery=stdev
      infs[i]->recov_time=int(round(norml(time_to_recovery, -dist_on_recovery, generator)));
  }

  infs[i]->quardt=100;infs[i]->testdt=100;//default no quarantining/testing
  if(randpercentage(quarp)){//to quarantine?
    if(dist_on_quardate>=0)
      infs[i]->quardt=(int)quardate+choosefrombin((int)dist_on_quardate);
    else
      infs[i]->quardt=int(round(norml(quardate, -dist_on_quardate, generator)));

    if(randpercentage(testp)){// to test?
      if(testdelay==0 || testdelay_shp<0)
	infs[i]->testdt=infs[i]->quardt + testdelay;//testing on fixed day after quarantine date
      else//testing delay follows a gamma distribution
	infs[i]->testdt=infs[i]->quardt+int(round(gamma(testdelay_shp, testdelay/testdelay_shp, generator)));
    }
  }

  //last operation (one greater than last operation)
  if(infs[i]->ill==-1){//dies (last op. is testing or death)
    infs[i]->lastop_time=infs[i]->dth_time;
    if(infs[i]->testdt!=100 && infs[i]->testdt > infs[i]->lastop_time)
      infs[i]->lastop_time=infs[i]->testdt;
  }
  else{//recovers (last op. is testing, death or seroconversion)
    infs[i]->lastop_time=infs[i]->recov_time;
    if(infs[i]->testdt!=100 && infs[i]->testdt > infs[i]->lastop_time)
      infs[i]->lastop_time=infs[i]->testdt;
    if(infs[i]->sero_time > infs[i]->lastop_time)
      infs[i]->lastop_time=infs[i]->sero_time;
  }
  (infs[i]->lastop_time)++;

  //set infection times
  if(inf_gam)//gamma distributed
    infs[i]->setinftimes(inf_tm_shp, inf_scl);
  else
    infs[i]->setinftimes(inf_start, inf_end);
  
  inflist[i]=1;
  return i;

}

void die(inf *a){//clear list position and delete
  inflist[a->num]=0;
  delete a;
  return;
}



int main(int argc, char *argv[]){
  int timeint;
  time_t timepoint;
  int i, tmpi, j, m, r, cur, num_runs;//number of runs
  int maxP;
  double R0;
  double trueR0,actualR0;
  int *P;
  int totdays;//total simulation length
  inf **infs=infar(0, MAXINFS-1);
  int init_infs;
  double avdthtime, avrecovtime, avtesttime, avserotime;
  int numdeaths, newdeaths, numrecovs;
  float dthrate;//percentage. A key parameter
  int geometric;//geometric or poisson or gamma? (geometric = 1, poisson = 0, gamma = -1)
  int gamswtch=0;
  // in case of gamma distribution
  double infscl=1.0;//scale for num to infect distribution
  double infshp;//shape for num to infect distribution
  int numinf;//number infected (cumulative)
  int numcurinf;//number currently infected
  int numcurinfold=0;
  int numinfectious;// number in the infectious window
  int newinfs; // number of new infections this time step. 
  int numquar;//number quarantined (cumulative)
  int numtest,newtests;//number tested (cumulative) and new
  int numill;//number ill (cumulative)
  int numsero;//cumulative seroconversion
  double percill;//percentage who fall (seriously) ill
  double percdeath;//percentage of ill who die
  //physical distancing?
  int haspd;//boolean
  int pd_at_dth;//pd starts at nth death
  int pd_at_test;//pd starts at nth tested infection
  int pd_at_inf;//pd starts at nth infection
  float pdeff1;//effectiveness of physical distancing
  float pdeff;//effectiveness of physical distancing. E.g. 40% - removes 2 in 5 contacts
  int pd;//physical distancing is occurring
  int inf_gam;//to gamma distribute infection times or not
  int inf_start, inf_end; //start and end of infective window
  double inf_mid, inf_tm_shp; // mean and shape parameter if gamma distributed
  double time_to_death, time_to_recovery, time_to_sero;//self explanatory
  double dist_on_death, dist_on_recovery, dist_on_sero;//binomial distributions: values 0,2,4,6
  double quarp;//percentage who get quarantined
  double testp;//percentage *of those quarantined* who are tested
  double quardate;//Currently assume all tests occur on a particular day in the infection cycle. Only those tested are quarantined. India: 10? UK: 12?
  double dist_on_quardate;//distribution on quardate
  double testdelay, testdelay_shp;

  double totpop;// total population (only relevant if herd=1)
  double effpop;//effective population (only relevant if herd=1)

  int haslockdown;//lockdown?
  int lockdownlen, lockdown2len;//length of lockdown
  int lockdownday, lockdown2day, lockdown2startday;//how many days into lockdown?
  float pdeff_lockdown, pdeff_lockdown2;//effectiveness of physical distancing post lockdown
  int lockdown_at_dth;//The lockdown begins after the death number lockdown_at_dth. UK ~200, India ~10
  int lockdown_at_test;//The lockdown begins after test number lockdown_at_test.
  int lockdown_at_inf;//The lockdown begins after infection number lockdown_at_inf.

  double popleak, popleak2;//leak into effective population post lockdown (an absolute value at the moment)
  int popleak_start_day=0, popleak2_start_day=0, popleak_end_day, popleak2_end_day;
  float infectible_proportion, infectible_proportion2;//default infectible proportion at lockdown

  int herd;//herd immunity?
  double herdlevel=0;
  char paramfilename[200], outfilename[200], endfname[206], logfname[204], datafilename[200];
  char tempword[200];
  FILE *fd, *fd1, *fd5, *fd6, *fd7; //files to store output
  int **alloutput;//to store all the simulation output
  double **avoutput, **SEoutput;//to store average, SE of output, synchronised
  double tmpSD;
  int *delays, totsims, maxdel;
  //FILE *fd3;

  //For the purposes of synchronising with data
  int syncflag=0;
  int sync_at_test;//at test number
  int sync_at_death;//at death number
  int sync_at_inf;//at infection number
  int sync_at_time;
  int syncclock=0;//start counting after synchronisation

  //These parameters are relevant if we want to 
  //run simulations upto or a certain number of days
  //beyond a particular death trigger
  //Currently hard-wired in, to avoid overloading parameter files
  int topresent=0;//only simulate to a fixed day namely "presentday" days after "trigger_dths" deaths or "trigger_infs" infections
  int trigger_dths=1;//Number of deaths which trigger the clock. Not for synchronisation
  int trigger_infs=0;//Number of infections which trigger the clock. Not for synchronisation
  int startinfs, endinfs;
  int presentday=1;//Number of days to run after trigger
  int startclock=0;//The clock

  //for doubling times
  int totdoubling=0;
  double avdoubling=0;
  double avinfs, avdths;//average infections and deaths at trigger point

  //
  int multiplier=1;
  int dynmultiply;//dynamic to speed up computation
  int cur_exp=1;
  int scale_at_infs;

  //data file
  int maxdat=1000, totdata=0;
  int **realdata=imatrix(0, maxdat-1, 0, 2);


  if(argc < 2){
    fprintf(stderr, "ERROR: you must provide a parameter file name. You may also provide an output file name.\n");
    exit(0);
  }
  strncpy (paramfilename, argv[1], sizeof(paramfilename));
  if(argc>=3)
    strncpy (outfilename, argv[2], sizeof(outfilename));
  else
    strcpy(outfilename, "data1/tmp");//default output file

  fd1=openftowrite(outfilename); //tab separated output
  strcpy(endfname, outfilename);strcat(endfname, "_av");
  fd5=openftowrite(endfname); //tab separated output - average values
  strcpy(endfname, outfilename);strcat(endfname, "_sync1");
  fd6=openftowrite(endfname); //tab separated output - average values
  strcpy(endfname, outfilename);strcat(endfname, "_sync");
  fd7=openftowrite(endfname); //tab separated output - values after synchronisation

  strcpy(logfname, outfilename);strcat(logfname, "_log");
  fd=openftowrite(logfname); //log file

  //options: general
  num_runs=getoptioni(paramfilename, "number_of_runs", 10, fd1);//model runs
  dthrate=getoptionf(paramfilename, "death_rate", 0.5, fd1);//death rate
  geometric=getoptioni(paramfilename, "geometric", 0, fd1);//default is Poisson distribution
  R0=getoptionf(paramfilename, "R0", 3.5, fd1);//basic reproduction number (approximately)
  infshp=getoptionf(paramfilename, "infshp", 0.1, fd1);//shape param
  totdays=getoptioni(paramfilename, "totdays", 150, fd1);//total simulation length
  totpop=getoptionf(paramfilename, "population", 66000000, fd1);//population
  inf_gam=getoptioni(paramfilename, "inf_gam", 0, fd1);//use gamma distribution for infection times? Default is no
  inf_start=getoptioni(paramfilename, "inf_start", 2, fd1);//start of infective window
  inf_end=getoptioni(paramfilename, "inf_end", 9, fd1);//end of infective window
  // if infection times are gamma distributed
  inf_mid=getoptionf(paramfilename, "inf_mid", 6, fd1);//mean infection time
  inf_tm_shp=getoptionf(paramfilename, "inf_tm_shp", 4, fd1);//shape parameter for infection time

  time_to_death=getoptionf(paramfilename, "time_to_death", 17, fd1);//survival time
  dist_on_death=getoptionf(paramfilename, "dist_on_death", -3, fd1);//distribution on time_to_death. Default = none
  time_to_recovery=getoptionf(paramfilename, "time_to_recovery", 20, fd1);//recovery time
  dist_on_recovery=getoptionf(paramfilename, "dist_on_recovery", -2, fd1);//distribution on time_to_recovery
  time_to_sero=getoptionf(paramfilename, "time_to_sero", 14, fd1);//seroconversion time
  dist_on_sero=getoptionf(paramfilename, "dist_on_sero", -3, fd1);//distribution on time_to_sero

  init_infs=getoptioni(paramfilename, "initial_infections", 10, fd1);//initial number infected
  herd=getoptioni(paramfilename, "herd", 1, fd1);//herd immunity?
  //options: quarantine and testing
  quarp=getoptionf(paramfilename, "percentage_quarantined", 4, fd1);//percentage of infecteds who are quarantined
  testp=getoptionf(paramfilename, "percentage_tested", 100, fd1);//the percentage *of those quarantined* who are tested
  if(getoption(paramfilename, "quardate", 1, tempword, 200)==0)
    quardate=getoptionf(paramfilename, "quardate", 12, fd1);//mean date of testing and quarantining
  else//legacy
    quardate=getoptionf(paramfilename, "testdate", 12, fd1);//mean date of testing and quarantining
  if(getoption(paramfilename, "dist_on_quardate", 1, tempword, 200)==0)
    dist_on_quardate=getoptionf(paramfilename, "dist_on_quardate", -3, fd1);//distribution on quarantine date
  else//legacy
    dist_on_quardate=getoptionf(paramfilename, "dist_on_testdate", -3, fd1);//distribution on quarantine date
  testdelay=getoptionf(paramfilename, "testdelay", 0, fd1);//mean delay from quarantining to testing
  testdelay_shp=getoptionf(paramfilename, "testdelay_shp", -1, fd1);//distribution on delay between quarantining and testing
  //options: lockdown
  haslockdown=getoptioni(paramfilename, "haslockdown", 0, fd1);//lockdown?
  lockdown_at_dth=getoptioni(paramfilename, "lockdth", -1, fd1);//lockdown at nth death
  lockdown_at_test=getoptioni(paramfilename, "lockdown_at_test", -1, fd1);//lockdown at nth test
  lockdown_at_inf=getoptioni(paramfilename, "lockdown_at_inf", -1, fd1);//lockdown at nth test
  lockdownlen=getoptioni(paramfilename, "lockdownlen", 0, fd1);//length of lockdown
  infectible_proportion=getoptionf(paramfilename, "infectible_proportion", 0.05555, fd1);
  pdeff_lockdown=getoptionf(paramfilename, "pdeff_lockdown", 60, fd1);//effectiveness of physical distancing after lockdown
  popleak=getoptionf(paramfilename, "popleak", 0, fd1);//leak into infectible population per day
  popleak_start_day=getoptioni(paramfilename, "popleak_start_day", 0, fd1);//when does the infectible population start to grow? The nth day of lockdown
  popleak_end_day=getoptioni(paramfilename, "popleak_end_day", 1000, fd1);//when does the infectible population end growing? Default is never.

  if(haslockdown==2){
    lockdown2startday=getoptioni(paramfilename, "lockdown2startday", 0, fd1);//start day of second lockdown
    lockdown2len=getoption2i(paramfilename, "lockdownlen", 0, fd1);//length of lockdown
    infectible_proportion2=getoption2f(paramfilename, "infectible_proportion", 0.05555, fd1);
    pdeff_lockdown2=getoption2f(paramfilename, "pdeff_lockdown", 60, fd1);//effectiveness of physical distancing after lockdown
    popleak2=getoption2f(paramfilename, "popleak", 0, fd1);//leak into infectible population per day
    popleak2_start_day=getoption2i(paramfilename, "popleak_start_day", 0, fd1);//when does the infectible population start to grow? The nth day of lockdown
    popleak2_end_day=getoption2i(paramfilename, "popleak_end_day", 1000, fd1);//when does the infectible population end growing? Default is never.
  }

  //options: physical distancing
  haspd=getoptioni(paramfilename, "physical_distancing", 0, fd1);//physical distancing?
  pd_at_dth=getoptioni(paramfilename, "pddth", -1,fd1);//physical distancing at nth death
  pd_at_test=getoptioni(paramfilename, "pd_at_test", -1,fd1);//physical distancing at nth recorded infection
  pd_at_inf=getoptioni(paramfilename, "pd_at_inf", -1,fd1);//physical distancing at nth infection
  pdeff1=getoptionf(paramfilename, "pdeff1", 30, fd1);//effectiveness of physical distancing

  //synchronisation with data
  sync_at_test=getoptionf(paramfilename, "sync_at_test", -1, fd1);//for synchronisation
  sync_at_inf=getoptionf(paramfilename, "sync_at_inf", -1, fd1);//for synchronisation
  sync_at_death=getoptionf(paramfilename, "sync_at_death", -1, fd1);//for synchronisation
  sync_at_time=getoptionf(paramfilename, "sync_at_time", -1, fd1);//for synchronisation
  // dynamic speeding up. Set to -1 for no speeding up
  scale_at_infs=getoptioni(paramfilename, "scale_at_infs", 50000,fd1);//default is to begin scaling at the 50000th infection

  if (getoption(paramfilename, "datafile", 1, datafilename, 200)==0){
    totdata=readDataFile(datafilename, realdata, maxdat);
  }

  if(scale_at_infs>0)
    dynmultiply=1;
  else
    dynmultiply=0;

  alloutput=imatrix(0, num_runs*totdays-1, 0, 9);
  avoutput=dmatrix(0, totdays-1, 0, 12);
  SEoutput=dmatrix(0, totdays-1, 0, 12);
  delays=(int *)malloc((size_t) ((num_runs)*sizeof(int)));

  //gamma distribution on individual R0 values
  if(geometric==-1){gamswtch=1;}
  infscl=R0/infshp;

  //fd3=fopen("data1/graph", "w");

  //How is the number to infect distributed? How to truncate?
  if(geometric==1){maxP=2*R0*(R0+1)<MAXDISCPROB?(int)(2*R0*(R0+1)):MAXDISCPROB;P=discGeom(R0, maxP);}//geometric (2 SD)
  else{maxP=3*R0<MAXDISCPROB?(int)(3*R0):MAXDISCPROB;P=discPois(R0, maxP);}//Poisson (three SD)

  percill=20.0;//percentage of people who fall quite ill (not currently used - for hospitalisations data?)
  percdeath=dthrate*100.0/percill;

  //random seeding
  timeint = time(&timepoint); /*convert time to an integer */
  srand(timeint);
  generator.seed(timeint);//seeding for distribution generator

  trueR0=0;
  if(!gamswtch){
    for(i=0;i<=maxP;i++){
      fprintf(fd, "%d\n", P[i]);
      if(i>=1){//true R0
	trueR0+=(i-1)*((double)(P[i-1]-P[i]))/1000.0;
      }
    }
    trueR0+=maxP*((double)(P[maxP]))/1000.0;
  }
  else{
    for(i=0;i<1000000;i++)
      trueR0+=round(gamma(infshp, infscl, generator))/1000000.0;
  }
  fprintf(fd, "R0=%.4f, trueR0=%.4f\n", R0, trueR0);
  fprintf(fd, "run\tactualR0\tavdthtime\tavrecovtime\tavtesttime\tavserotime\n");

  //nest order: For each run... for each day... for each individual
  avinfs=0.0;avdths=0.0;
  for(r=0;r<num_runs;r++){//Each model run
    startclock=0;
    numinf=0;numcurinf=0;numcurinfold=0;numdeaths=0;newdeaths=0;numrecovs=0;
    numquar=0;numtest=0;newtests=0;numill=0;numsero=0;
    actualR0=0;avdthtime=0;avrecovtime=0;avserotime=0;avtesttime=0;
    lockdownday=0;lockdown2day=0;
    effpop=totpop;
    pd=0;
    syncflag=0;syncclock=0;
    cur_exp=1;
    multiplier=1;

    for(i=0;i<MAXINFS;i++){inflist[i]=0;}//initialise

    for(i=0;i<init_infs;i++){
      cur=create(infs, gamswtch, infshp, infscl, P, maxP, inf_gam, inf_start, inf_end, inf_mid, inf_tm_shp, &numinf, &numcurinf, &newinfs, &numill, percill, percdeath, time_to_death, dist_on_death, time_to_recovery, dist_on_recovery, time_to_sero, dist_on_sero, quardate, quarp, dist_on_quardate, testp, testdelay, testdelay_shp);//create
      //fprintf(fd3, "0 %d\n", cur);

      actualR0=actualR0*((double)(numinf-1))/((double)(numinf))+(double)((infs[cur])->numtoinf)/((double)(numinf));

      //fprintf(stderr, "actualR0=%.4f\n", actualR0);
      fprintf(stderr, "infs[%d] (illstate=%d) will infect %d at times:\n",cur, infs[cur]->ill, infs[cur]->numtoinf);
      for(j=0;j<(infs[cur])->numtoinf;j++){
	fprintf(stderr, "   %d\n", infs[cur]->inftimes[j]);
      }
    }
    
    for(m=0;m<totdays;m++){//each day

      //rescale. kill off half randomly; double weight of remainder
      if(dynmultiply && numcurinf>=scale_at_infs*int_pow(2,cur_exp-1) && multiplier==int_pow(2,cur_exp-1)){//multiplier
	cur_exp++;multiplier*=2;
	for(i=0;i<MAXINFS;i++){// kill off every other active infection
	  if(inflist[i]==1 && randnum(2)<1)
	    die(infs[i]);//no removal from stats
	}
      }

      numinfectious=0;newinfs=0;newdeaths=0;newtests=0;
      numcurinfold=numcurinf;
      if(herd)
	fprintf(stderr, "herdlevel=%.4f\n", herdlevel);

      if(haspd && ((pd_at_dth>0 && numdeaths>=pd_at_dth) || (pd_at_test>0 && numtest>=pd_at_test) || (pd_at_inf>0 && numinf>=pd_at_inf))){//physical distancing
	pd=1;
	pdeff=pdeff1;
      }
      else
	pd=0;

      if(haslockdown){
	//lockdown 2 takes precedence
	if(haslockdown==2 && lockdownday>=lockdown2startday && lockdown2day<lockdown2len){
	  if(lockdown2day==0){
	    effpop=totpop;
	    effpop*=infectible_proportion2;//effective infectible population drops
	    fprintf(stderr, "Lockdown 2 starts. Effective population now %.4f.\n", effpop);
	  }
	  else{ 
	    if(lockdown2day>=popleak2_start_day && lockdown2day<=popleak2_end_day)
	      effpop+=popleak2;//leak into infectible population
	    fprintf(stderr, "In lockdown 2. Effective population now %.4f.\n", effpop);
	  }
	  pd=1;
	  pdeff=pdeff_lockdown2;//physical distancing becomes more effective  
	  lockdown2day++;
	}
	else if(((lockdown_at_dth>0 && numdeaths>=lockdown_at_dth) || (lockdown_at_test>0 && numtest>=lockdown_at_test) || (lockdown_at_inf>0 && numinf>=lockdown_at_inf)) && lockdownday<lockdownlen){
	  if(lockdownday==0){
	    effpop*=infectible_proportion;//effective infectible population drops
	    fprintf(stderr, "Lockdown 1 starts. Effective population now %.4f.\n", effpop);
	  }
	  else{ 
	    if(lockdownday>=popleak_start_day && lockdownday<=popleak_end_day)
	      effpop+=popleak;//leak into infectible population
	    fprintf(stderr, "In lockdown 1. Effective population now %.4f.\n", effpop);
	  }
	  pd=1;
	  pdeff=pdeff_lockdown;//physical distancing becomes more effective  
	  lockdownday++;
	}
	else{//lockdown finishes. Assume physical distancing returns to early levels
	  effpop=totpop;
	  if(lockdownday>=lockdownlen){
	    fprintf(stderr, "Lockdown finished. Effective population now %.4f.\n", effpop);
	    lockdownday++;//to know when to enter lockdown2
	  }
	  if(haspd){
	    pdeff=pdeff1;
	  }
	  
	}
      }
      if(pd){
	fprintf(stderr, "physical distancing = %.2f.\n", pdeff);
      }


      for(i=0;i<MAXINFS;i++){//for each infected person
	if(inflist[i]==1){
	  (infs[i]->age)++;//age updates at start...
	  if(infs[i]->age==infs[i]->lastop_time){//done with
	    die(infs[i]);//deallocate
	    continue;
	  }

	  if(infs[i]->age >= inf_start && infs[i]->age <= inf_end){//so far kept as is regardless of distribution
	    numinfectious+=multiplier;
	  }

	  if(infs[i]->age==infs[i]->sero_time){//seroconversion
	    numsero+=multiplier;
	    avserotime=avserotime*((double)(numsero-multiplier))/((double)(numsero))+(double)(multiplier*(infs[i])->age)/((double)(numsero));
	  }

	  if(infs[i]->age==infs[i]->quardt){//quarantine?
	    infs[i]->quar=1;
	    numquar+=multiplier;
	  }

	  if(infs[i]->age==infs[i]->testdt){//test?
	    numtest+=multiplier;newtests+=multiplier;
	    avtesttime=avtesttime*((double)(numtest-multiplier))/((double)(numtest))+(double)(multiplier*(infs[i])->age)/((double)(numtest));
	  }

	  if(infs[i]->ill==-1 && infs[i]->age==infs[i]->dth_time){//die
	    numdeaths+=multiplier;newdeaths+=multiplier;numcurinf-=multiplier;
	    avdthtime=avdthtime*((double)(numdeaths-multiplier))/((double)(numdeaths))+(double)(multiplier*(infs[i])->age)/((double)(numdeaths));
	  }
	  else if(infs[i]->ill!=-1 && infs[i]->age==infs[i]->recov_time){//recover
	    numcurinf-=multiplier;numrecovs+=multiplier;
	    avrecovtime=avrecovtime*((double)(numrecovs-multiplier))/((double)(numrecovs))+(double)(multiplier*(infs[i])->age)/((double)(numrecovs));
	  }
	  else if(infs[i]->quar==0 && infs[i]->age<MAXAGE && (!pd|| (pd && randpercentage(100.0-pdeff)))){//still being processed, not quarantined, no physical distancing or pd not happening
	    if(herd){herdlevel=100.0*((double)numinf/(double)effpop);}
	    if(!herd || (herd && randpercentage(100.0-herdlevel))){
	      for(j=0;j<infs[i]->infnums[infs[i]->age];j++){
		tmpi=create(infs, gamswtch, infshp, infscl, P, maxP, inf_gam, inf_start, inf_end, inf_mid, inf_tm_shp, &numinf, &numcurinf, &newinfs, &numill, percill, percdeath, time_to_death, dist_on_death, time_to_recovery, dist_on_recovery, time_to_sero, dist_on_sero, quardate, quarp, dist_on_quardate, testp, testdelay, testdelay_shp);//create new infecteds
		numinf+=(multiplier-1);numcurinf+=(multiplier-1);newinfs+=(multiplier-1);
		actualR0=actualR0*((double)(numinf-multiplier))/((double)(numinf))+(double)(multiplier*(infs[tmpi])->numtoinf)/((double)(numinf));
		if(infs[tmpi]->ill==1)
		  numill+=(multiplier-1);
		//		    fprintf(fd3, "%d %d %d\n%d %d %d\n\n", m, i, i, m+1, tmpi, tmpi);
		if(tmpi>i)//this will update to zero as they come later in the sequence
		  infs[tmpi]->age--;
	      }
	    }
	  }
	}
      }//cycled through all infected individuals

      fprintf(stderr, "%d: numinf=%d, newinfs=%d, numcurinf=%d(%.2fpc), numdeaths=%d, newdeaths=%d, numtest=%d, numinfectious=%d, numsero=%d\n", m, numinf, newinfs, numcurinf, numcurinfold>=1?100.0*((double)numcurinf-(double)numcurinfold)/((double)numcurinfold):-1,numdeaths, newdeaths, numtest, numinfectious, numsero);
      fprintf(fd1,"%d\t%d\t%d\t%d\t %d\t%d\t%d\t%d\t%d\t%d\n", m, numinf, newinfs, numcurinf, numdeaths, newdeaths, numtest,newtests,numinfectious,numsero);
      alloutput[r*totdays+m][0]=m;alloutput[r*totdays+m][1]=numinf;
      alloutput[r*totdays+m][2]=newinfs;alloutput[r*totdays+m][3]=numcurinf;
      alloutput[r*totdays+m][4]=numdeaths;alloutput[r*totdays+m][5]=newdeaths;
      alloutput[r*totdays+m][6]=numtest;alloutput[r*totdays+m][7]=newtests;
      alloutput[r*totdays+m][8]=numinfectious;alloutput[r*totdays+m][9]=numsero;

      //Setting the delays
      if(!syncflag && ((sync_at_test>0 && numtest>=sync_at_test) || (sync_at_death>0 && numdeaths>=sync_at_death) || (sync_at_inf>0 && numinf>=sync_at_inf))){
	delays[r]=m-sync_at_time;
	syncflag=1;
      }

      //      fprintf(fd3, "\n");fflush(fd3);
      fflush(fd);fflush(fd1);

      //Are we running the model only to a particular moment?
      if(topresent){
	if(trigger_infs){//triggered by infection numbers
	  if(numinf>=trigger_infs){
	    //printf("%d, %d: numinf=%d, numdeaths=%d, \n", r, m, numinf, numdeaths);
	    if(startclock==0)
	      startinfs=numinf;
	    startclock++;
	  }
	  if(startclock==presentday+1){
	    //printf("%d, %d: numinf=%d, numdeaths=%d, \n", r, m, numinf, numdeaths);
	    endinfs=numinf;
	    printf("model run %d: %d %d %d\n", r, startinfs, endinfs, presentday);
	    if(presentday>0 && endinfs-startinfs!=0){
	      printf("doubling=%.4f\n", log(2.0)*(presentday)/(log(endinfs)-log(startinfs)));
	      totdoubling++; avdoubling+=log(2.0)*(presentday)/(log(endinfs)-log(startinfs));
	    }
	    avinfs+=numinf;avdths+=numdeaths;
	    break;
	  }
	}
	else{
	  if(numdeaths>=trigger_dths){
	    if(startclock==0)
	      startinfs=numinf;
	    startclock++;
	  }
	  
	  if(startclock==presentday+1){
	    endinfs=numinf;
	    printf("model run %d: %d %d %d\n", r, startinfs, endinfs, presentday);
	    if(presentday>0 && endinfs-startinfs!=0){
	      printf("doubling=%.4f\n", log(2.0)*(presentday)/(log(endinfs)-log(startinfs)));
	      totdoubling++; avdoubling+=log(2.0)*(presentday)/(log(endinfs)-log(startinfs));
	    }
	    avinfs+=numinf;avdths+=numdeaths;
	    break;
	  }
	}
      }
      if(syncflag)
	syncclock++;// counts days since synchronisation event
    }
    if(!syncflag){//synchronisation point never reached (died out?)
      delays[r]=0;
      syncflag=1;
    }

    //Only output to synchronisation file if there is a data file and synchronisation point reached and positive delay
    if(delays[r]>0){
      for(m=0;m<totdays-delays[r];m++){
	for(i=0;i<10;i++)
	  fprintf(fd7, "%d\t", alloutput[r*totdays+m+delays[r]][i]);
	if(m<totdata){
	  for(i=0;i<3;i++)
	    fprintf(fd7, "%d\t", realdata[m][i]);
	}
	else{
	  for(i=0;i<3;i++)
	    fprintf(fd7, "?\t");
	}
	fprintf(fd7, "\n");
      }
      fprintf(fd7, "\n");
      fflush(fd7);
    }

    fprintf(fd1,"\n");
    for(i=0;i<MAXINFS;i++){//free memory
      if(inflist[i]==1)
	die(infs[i]);//deallocate (numcurinf will get reset anyway)
    }
    fprintf(fd, "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", r+1, actualR0, avdthtime, avrecovtime, avtesttime, avserotime);
  }

  totsims=0;
  maxdel=0;
  for(r=0;r<num_runs;r++){
    if(delays[r]>0){
      totsims++;
      if(delays[r]>maxdel){maxdel=delays[r];}
    }
  }
  //only output to average file if there is a data file and synchronisation point reached and positive delay

  for(m=0;m<totdays-maxdel;m++){
    for(i=0;i<10;i++){//average values
      avoutput[m][i]=0;
      for(r=0;r<num_runs;r++){
	if((totsims>0 && delays[r]>0) || totsims==0)
	  avoutput[m][i]+=alloutput[r*totdays+m+delays[r]][i];
      }
      if(totsims>0)
	avoutput[m][i]/=totsims;
      else
	avoutput[m][i]/=num_runs;
      fprintf(fd5, "%.1f\t", avoutput[m][i]);
     
    }
    for(i=0;i<10;i++){//standard errors
      tmpSD=0.0;
      for(r=0;r<num_runs;r++){
	if((totsims>0 && delays[r]>0) || totsims==0)
	  tmpSD+=(alloutput[r*totdays+m+delays[r]][i]-avoutput[m][i])*(alloutput[r*totdays+m+delays[r]][i]-avoutput[m][i]);
      }
      if(totsims>1){
	tmpSD/=((double)totsims-1.0);//population SD
	SEoutput[m][i]=sqrt(tmpSD/(double)totsims);//SE
      }      
      else if(num_runs>1){
	tmpSD/=((double)num_runs-1.0);
	SEoutput[m][i]=sqrt(tmpSD/(double)num_runs);
      }
      else{
	SEoutput[m][i]=0.0;
      }
      fprintf(fd5, "%.4f\t", SEoutput[m][i]);
    }

    if(m<totdata){//output data from the datafile too
      for(i=0;i<3;i++)
	fprintf(fd5, "%d\t", realdata[m][i]);
    }
    else{
      for(i=0;i<3;i++)
	fprintf(fd5, "?\t");
    }
    fprintf(fd5, "\n");

    for(r=0;r<num_runs;r++){//grouped sync file
      fprintf(fd6, "%d\t", m);
      for(i=1;i<10;i++){
	if((totsims>0 && delays[r]>0) || totsims==0)
	  fprintf(fd6, "%d\t", alloutput[r*totdays+m+delays[r]][i]);
	else
	  fprintf(fd6, "0\t");
      }
      if(m<totdata){//output data from the datafile too
	for(i=0;i<3;i++)
	  fprintf(fd6, "%d\t", realdata[m][i]);
      }
      else{
	for(i=0;i<3;i++)
	  fprintf(fd6, "?\t");
      }
      fprintf(fd6, "\n");
    }
    fprintf(fd6, "%d\t", m);
    for(i=1;i<10;i++)
      fprintf(fd6, "%.1f\t", avoutput[m][i]);

    if(m<totdata){//death undercount, mean + 95%CI
      fprintf(fd6, "%.4f\t", 100.0*(avoutput[m][4]-realdata[m][1])/avoutput[m][4]);
      fprintf(fd6, "%.4f\t", 100.0*((avoutput[m][4]-1.96*SEoutput[m][4])-realdata[m][1])/(avoutput[m][4]-1.96*SEoutput[m][4]));
      fprintf(fd6, "%.4f\t", 100.0*((avoutput[m][4]+1.96*SEoutput[m][4])-realdata[m][1])/(avoutput[m][4]+1.96*SEoutput[m][4]));
    }

    fprintf(fd6, "\n");

    fprintf(fd6, "%d\t", m);
    for(i=1;i<10;i++)
      fprintf(fd6, "%.4f\t", SEoutput[m][i]);

    if(m<totdata){//case undercount, mean + 95%CI
      fprintf(fd6, "%.4f\t", 100.0*(avoutput[m][6]-realdata[m][0])/avoutput[m][6]);
      fprintf(fd6, "%.4f\t", 100.0*((avoutput[m][6]-1.96*SEoutput[m][6])-realdata[m][0])/(avoutput[m][6]-1.96*SEoutput[m][6]));
      fprintf(fd6, "%.4f\t", 100.0*((avoutput[m][6]+1.96*SEoutput[m][6])-realdata[m][0])/(avoutput[m][6]+1.96*SEoutput[m][6]));
    }
    fprintf(fd6, "\n\n");

  }

  for(m=totdays-maxdel;m<totdays;m++){
    for(i=0;i<10;i++){
      avoutput[m][i]=-1;
      fprintf(fd5, "%.1f\t", avoutput[m][i]);
    }
    for(i=0;i<10;i++){
      SEoutput[m][i]=0;
      fprintf(fd5, "%.1f\t", SEoutput[m][i]);
    }
    for(i=0;i<3;i++)
      fprintf(fd5, "?\t");
    fprintf(fd5, "\n");
  }

  if(topresent){
    if(totdoubling>0)
      printf("avinfs=%.4f, avdeaths=%.4f, av. doubling time=%.4f\n", avinfs/((double)num_runs), avdths/((double)num_runs),avdoubling/((double)totdoubling));
    else
      printf("avinfs=%.4f, avdeaths=%.4f\n", avinfs/((double)num_runs), avdths/((double)num_runs));
  }


  free_infar(infs, 0, MAXINFS-1);
  free((char *) P);//fclose(fd3);
  fclose(fd);fclose(fd1);fclose(fd5);fclose(fd6);fclose(fd7);
  free_imatrix(alloutput, 0, num_runs*totdays-1, 0, 9);
  free_dmatrix(avoutput, 0, totdays-1, 0, 9);
  free_dmatrix(SEoutput, 0, totdays-1, 0, 9);
  free_imatrix(realdata, 0, maxdat-1, 0, 2);
  free((char*)delays);
  return 0;
}

