/* Copyright (C) 2020, Murad Banaji
 *
 * This file is part of COVIDSIM
 *
 * COVIDSIM is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 3, 
 * or (at your option) any later version.
 *
 * COVIDSIM is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with COVIDSIM: see the file COPYING.  If not, see 
 * <https://www.gnu.org/licenses/>

 */


#include "inf.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> // random seeding
#include <ctype.h>

// Maximum number of infected individuals - memory limit?
#define MAXINFS 10000000

//Externally declared (bad practice I know!)

int inflist[MAXINFS]; //For book-keeping free spaces in list

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

int getnthblock(char *s, char *v, int len, int n){
  // get the nth valid block (only spaces count as separators) from a string s and put it in v. Returns the next position in  s. 
  int i, j, k;
  i=0, k=0;
  if(len < 2){
    fprintf(stderr, "ERROR in getnthblock in stringmanip.c: third argument must be at least 2.\n");
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
      fprintf(stderr, "WARNING: in routine getnthblock in file stringmanip.c: word is longer than maximum length.\n");
    }
    return k;
  }
  return 0;
}


int getoption(char *fname, const char optname[], char v[], int max){
  FILE *fd;
  int len, j, flag=0;
  char oneline[200];
  char modname[50];
  int lim=200;
  fd = fopen(fname, "r");
  if(fd == NULL){
    fprintf(stderr, "ERROR in routine getoption in \"stringmanip.c\": file %s could not be opened for reading.\n", fname);
    return -1;
  }
  while((len = getline(fd, oneline, lim)) > 0){
    j=0;
    while((isspace((int) oneline[j])) || (oneline[j] == 13)){j++;}
    if ((oneline[j] == '/') || (oneline[j] == '\n') || (oneline[j] == '\0')){} // comment/empty lines
    else{
      getnthblock(oneline, modname, 50, 1);
      if(strcmp(modname, optname) == 0){
        flag = 1;
        getnthblock(oneline, modname, 50, 2);
        if((int)(strlen(modname)) < max-1)
	  strcpy(v, modname);
        else{
	  fprintf(stderr, "ERROR in routine getoption in \"stringmanip.c\": Option %s in file %s has value %s which is too long.\n", optname, fname, modname);
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
    fprintf(stderr, "ERROR in routine getoption in \"stringmanip.c\": Option %s could not be found in file %s.\n", optname, fname);
    v[0] = '\0';
    return -2;
  }
  return 0;

}

int getoptioni(char *fname, const char optname[], int defval, FILE *fd1){
  char tempword[200];
  int val;
  if(getoption(fname, optname, tempword, 200)!=0)
    val=defval;
  else
    val=atoi(tempword);
  fprintf(fd1, "#%s %d\n", optname, val);
  return val;
}

float getoptionf(char *fname, const char optname[], float defval, FILE *fd1){
  char tempword[200];
  float val;
  if(getoption(fname, optname, tempword, 200)!=0)
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

//Setting the tail to zero in the distributions below lowers the true expected value. Also given the different weights of the tails, it can make the dynamics for the two distributions differ.

//Poisson distribution cast to integers (not rounded) and with the tail given to zero 
int *discPoisold(double lamb, int max){
  int *P;
  int k, Ptot;
  P=(int *)malloc((size_t) ((max+1)*sizeof(int)));
  Ptot=0;
  for(k=max;k>0;k--){
    Ptot+=(int)(1000.0*Pois(lamb, k));
    P[k]=Ptot;
  }
  P[0]=1000;
  return P;
}

//Poisson distribution cast to integers (not rounded) and with the tail given to largest value
int *discPois(double lamb, int max){
  int *P;
  int k, Ptot=0;
  P=(int *)malloc((size_t) ((max+1)*sizeof(int)));
  P[0]=1000;
  for(k=1;k<=max;k++){
    P[k]=P[k-1]-(int)(1000.0*Pois(lamb, k-1));
  }
  return P;
}



//Geometric distribution cast to integers (not rounded) and with the tail given to zero
int *discGeomold(double lamb, int max){
  int *P;
  int k, Ptot;
  P=(int *)malloc((size_t) ((max+1)*sizeof(int)));
  Ptot=0;
  for(k=max;k>0;k--){
    Ptot+=(int)(1000.0*Geom(lamb, k));
    P[k]=Ptot;
  }
  P[0]=1000;
  return P;
}

//Geometric distribution cast to integers (not rounded) and with the tail given to zero
int *discGeom(double lamb, int max){
  int *P;
  int k, Ptot=0;
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

// infected class
// inf P is the distribution, maxP is the greatest num to infect
// time steps in days
// infected at age 0
// quarantined? Then can't infect
// ill? Can die if ill
// How many will this person infect? Follows the distribution
// The assumption is that it isn't just randomly about contact
// perhaps also biology?
// The infection times are chosen from a uniform distribution on some
// range of values


inf::inf(int orgnum, int P[], int maxP){
  age = 0;
  ill = 0;
  quar = 0;
  num = orgnum;
  numtoinf=choosefromdist(P, maxP);
}


// Set the times at which infection occurs
void inf::setinftimes(int rmin, int rmax){
  int i; 
  for(i=0;i<numtoinf;i++){
    inftimes[i]=randnum(rmax-rmin+1)+rmin;
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
/* This function allocates a set of character pointers */
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
/* free a float chararry allocated by charar() */
{
  free((char *) (m+nl-1));
}



int create(inf *infs[], int P[], int maxP, int inf_start, int inf_end, int *numinf, int *numcurinf, int *numill, double percill){
  int i=firstfreepos(inflist, MAXINFS);
  int j;
  if(i==-1)//no more space
    exit(0);
  infs[i] = new inf(i, P, maxP);
  (*numinf)++;
  (*numcurinf)++;

  for(j=0;j<21;j++){//number to infect at time j
    infs[i]->infnums[j]=0;
  }
  // who falls ill?
  if(randpercentage(percill)){
    infs[i]->ill=1;
    (*numill)++;
    //fprintf(stderr, "ill=%d\n", infs[i]->ill);
  }
  //set infection times
  infs[i]->setinftimes(inf_start, inf_end);
  inflist[i]=1;
  return i;

}

void die(inf *a, int *numcurinf){
  (*numcurinf)--;
  inflist[a->num]=0;
  delete a;
  return;
}


int main(int argc, char *argv[]){
  int timeint;
  time_t timepoint;
  int i, j, m, r, cur, num_runs;//number of runs
  int maxP;
  double R0;
  double trueR0;
  int *P;
  int totdays;//total simulation length
  inf **infs=infar(0, MAXINFS-1);
  int init_infs;
  int numdeaths;
  float dthrate;//percentage. A key parameter
  int geometric;//geometric or poisson?
  int numinf;//number infected (cumulative)
  int numcurinf;//number currently infected
  int numcurinfold=0;
  int numquar;//number quarantined (cumulative)
  int numtest;//number tested (cumulative)
  int numill;//number ill (cumulative)
  double percill;//percentage who fall (seriously) ill
  double percdeath;//percentage of ill who die
  int haslockdown;//lockdown?
  int lockdownlen;//length of lockdown
  int lockdownday;//how many days into lockdown?
  //physical distancing?
  int haspd;//boolean
  int pd_at_dth;//starts at nth death
  float pdeff1;//effectiveness of physical distancing
  float pdeff_lockdown;//effectiveness of physical distancing post lockdown
  float pdeff;//effectiveness of physical distancing. E.g. 40% - removes 2 in 5 contacts
  int pd;//physical distancing is occurring
  int startclock=0;//Used when we want to run a certain number of days after some trigger point
  double avinfs, avdths;
  int inf_start, inf_end; //start and end of infective window
  int time_to_death, time_to_recovery;//self explanatory
  int topresent=0;//only simulate to a fixed day, namely "presentday" days after "trigger_dths" deaths
  int trigger_dths=10;//Number of deaths which trigger the clock
  int presentday=18;//Number of days to run after trigger
  double quarp;//percentage who get quarantined
  double testp;//percentage *of those quarantined* who are tested
  int testdate;//Currently assume all tests occur on a particular day in the infection cycle. Only those tested are quarantined. India: 10? UK: 12?
  int lockdown_at_dth;//The lockdown begins after the death number lockdown_at_dth. UK ~200, India ~10
  double totpop;// total population (only relevant if herd=1)
  double effpop;//effective population (only relevant if herd=1)
  double popleak;//leak into effective population post lockdown (an absolute value at the moment)
  int popleak_start_day=0;
  float infectible_proportion;//default infected proportion at lockdown
  int herd;//herd immunity?
  double herdlevel=0;
  char paramfilename[200], outfilename[200], fname[200];
  FILE *fd, *fd1; //files to store output


  if(argc < 2){
    fprintf(stderr, "ERROR: you must provide a parameter file name. You may also provide an output file name.\n");
    exit(0);
  }
  strncpy (paramfilename, argv[1], sizeof(paramfilename));
  if(argc>=3)
    strncpy (outfilename, argv[2], sizeof(outfilename));
  else
    strcpy(outfilename, "data1/tmp");//default output file


  fd1=fopen(outfilename, "w"); //tab separated output

  //options: general
  num_runs=getoptioni(paramfilename, "number_of_runs", 10, fd1);//model runs
  dthrate=getoptionf(paramfilename, "death_rate", 0.5, fd1);//death rate
  geometric=getoptioni(paramfilename, "geometric", 0, fd1);//default is Poisson distribution
  R0=getoptionf(paramfilename, "R0", 3.5, fd1);//basic reproduction number (approximately)
  totdays=getoptioni(paramfilename, "totdays", 150, fd1);//total simulation length
  totpop=getoptionf(paramfilename, "population", 66000000, fd1);//population
  inf_start=getoptioni(paramfilename, "inf_start", 3, fd1);//start of infective window
  inf_end=getoptioni(paramfilename, "inf_end", 14, fd1);//end of infective window
  time_to_death=getoptioni(paramfilename, "time_to_death", 17, fd1);//survival time
  time_to_recovery=getoptioni(paramfilename, "time_to_recovery", 20, fd1);//recovery time
  init_infs=getoptioni(paramfilename, "initial_infections", 10, fd1);//initial number infected
  herd=getoptioni(paramfilename, "herd", 1, fd1);//herd immunity?
  //options: quarantine and testing
  quarp=getoptionf(paramfilename, "percentage_quarantined", 4, fd1);//percentage of infecteds who are quarantined
  testp=getoptionf(paramfilename, "percentage_tested", 100, fd1);//the percentage *of those quarantined* who are tested
  testdate=getoptioni(paramfilename, "testdate", 12, fd1);//leak into infectible population per day
  //options: lockdown
  haslockdown=getoptioni(paramfilename, "haslockdown", 0, fd1);//lockdown?
  lockdown_at_dth=getoptioni(paramfilename, "lockdth", 200, fd1);//lockdown at nth death
  lockdownlen=getoptioni(paramfilename, "lockdownlen", 0, fd1);//length of lockdown
  infectible_proportion=getoptionf(paramfilename, "infectible_proportion", 0.05555, fd1);
  pdeff_lockdown=getoptionf(paramfilename, "pdeff_lockdown", 60, fd1);//effectiveness of physical distancing after lockdown
  popleak=getoptionf(paramfilename, "popleak", 0, fd1);//leak into infectible population per day
  popleak_start_day=getoptioni(paramfilename, "popleak_start_day", 0, fd1);//when does the infectible population start to grow? The nth day of lockdown
  //options: physical distancing
  haspd=getoptioni(paramfilename, "physical_distancing", 0, fd1);//physical distancing?
  pd_at_dth=getoptioni(paramfilename, "pddth", 1,fd1);//physical distancing at nth death
  pdeff1=getoptionf(paramfilename, "pdeff1", 30, fd1);//effectiveness of physical distancing

  if(haslockdown)
    sprintf(fname, "data2/dth%.1f_%s_lamb%.2f_lock.csv", dthrate, geometric==1?"geom":"pois", R0);
  else
    sprintf(fname, "data2/dth%.1f_%s_lamb%.2f_nolock.csv", dthrate, geometric==1?"geom":"pois", R0);
  fprintf(stderr, "%s\n", fname);
  fd=fopen(fname, "w");//csv output

  //How is the number to infect distributed?
  if(geometric){maxP=1.5*R0*(R0+1)<30?(int)(1.5*R0*(R0+1)):30;P=discGeom(R0, maxP);}//geometric (1.5 SD)
  else{maxP=R0<10?(int)(3*R0):30;P=discPois(R0, maxP);}//Poisson (three SD)


  percill=20.0;//percentage of people who fall quite ill (not currently used - for hospitalisations data?)
  percdeath=dthrate*100.0/percill;

  //random seeding
  timeint = time(&timepoint); /*convert time to an integer */
  srand(timeint);

  trueR0=0;
  for(i=0;i<=maxP;i++){
    fprintf(stderr, "%d\n", P[i]);
    if(i>=1){//true R0
      trueR0+=(i-1)*((double)(P[i-1]-P[i]))/1000.0;
    }
  }
  trueR0+=maxP*((double)(P[maxP]))/1000.0;
  fprintf(stderr, "R0=%.4f, trueR0=%.4f\n", R0, trueR0);
  //exit(0);

  //At each time step go through each infected, augment its age, check its numtoinf, and then go through its inftimes, and spawn new if needed, remove it if needed.
  avinfs=0.0;avdths=0.0;
  for(r=0;r<num_runs;r++){//Each model run
    startclock=0;
    numinf=0;//cumulative total infections
    numcurinf=0;//current infections
    numcurinfold=0;
    numquar=0;
    numtest=0;
    numill=0;
    numdeaths=0;
    lockdownday=0;
    effpop=totpop;
    pd=0;
    for(i=0;i<MAXINFS;i++){inflist[i]=0;}//initialise

    for(i=0;i<init_infs;i++){
      cur=create(infs, P, maxP, inf_start, inf_end, &numinf, &numcurinf, &numill, percill);//create
      fprintf(stderr, "infs[%d] (illstate=%d) will infect %d at times:\n",cur, infs[cur]->ill, infs[cur]->numtoinf);
      for(j=0;j<(infs[cur])->numtoinf;j++){
	fprintf(stderr, "   %d\n", infs[cur]->inftimes[j]);
      }
    }
    
    for(m=0;m<totdays;m++){//each day
      numcurinfold=numcurinf;
      if(herd)
	fprintf(stderr, "herdlevel=%.4f\n", herdlevel);

      if(haspd && numdeaths>=pd_at_dth){//physical distancing
	pd=1;
	pdeff=pdeff1;
      }
      else
	pd=0;

      if(haslockdown){
	if(numdeaths>=lockdown_at_dth && lockdownday<lockdownlen){
	  if(lockdownday==0){
	    effpop*=infectible_proportion;//effective infectible population drops
	    fprintf(stderr, "Lockdown starts. Effective population now %.4f.\n", effpop);
	  }
	  else{ 
	    if(lockdownday>=popleak_start_day)
	      effpop+=popleak;//leak into infectible population
	    fprintf(stderr, "In lockdown. Effective population now %.4f.\n", effpop);
	  }
	  pd=1;
	  pdeff=pdeff_lockdown;//physical distancing becomes more effective  
	  lockdownday++;
	}
      }
      if(pd){
	fprintf(stderr, "physical distancing = %.2f.\n", pdeff);
      }


      for(i=0;i<MAXINFS;i++){//for each infected person
	if(inflist[i]==1){
	  (infs[i]->age)++;//age updates at start...?
	  if(infs[i]->age==testdate && randpercentage(quarp)){//quarantine?
	    infs[i]->quar=1;
	    numquar++;
	    if(randpercentage(testp))
	      numtest++;
	  }
	  if((infs[i]->age>time_to_death && infs[i]->ill==1) || infs[i]->age>time_to_recovery){
	    if((infs[i]->ill)==1){//if ill
	      if(randpercentage(percdeath))
		numdeaths++;
	    }
	    die(infs[i], &numcurinf);// die or recover
	  }
	  else if(infs[i]->quar==0){//if not quarantined, then infect
	    if(!pd|| (pd && randpercentage(100.0-pdeff))){//physical distancing
	      if(herd){herdlevel=100.0*((double)numinf/(double)effpop);}
	      if(!herd || (herd && randpercentage(100.0-herdlevel))){
		for(j=0;j<infs[i]->infnums[infs[i]->age];j++)
		  create(infs, P, maxP, inf_start, inf_end, &numinf, &numcurinf, &numill, percill);//create new infecteds
	      }
	    }
	  }
	}
      }//cycled through all infected individuals

      fprintf(stderr, "%d: numinf=%d, numcurinf=%d(%.2fpc), numdeaths=%d, numtest=%d, numill=%d\n", m, numinf, numcurinf, numcurinfold>=1?100.0*((double)numcurinf-(double)numcurinfold)/((double)numcurinfold):-1,numdeaths, numtest, numill);
      fprintf(fd,"%d,%d,%.4f,%d,%d,%.4f, %d,%.4f,%d\n", m, numinf, numinf>=1?log10((float)numinf):-1, numcurinf, numdeaths, numdeaths>=1?log10((float)numdeaths):-1,numtest,numtest>=1?log10((float)numtest):-1, numill);
      fprintf(fd1,"%d\t%d\t%.4f\t%d\t%d\t%.4f\t %d\t%.4f\t%d\n", m, numinf, numinf>=1?log10((float)numinf):-1, numcurinf, numdeaths, numdeaths>=1?log10((float)numdeaths):-1,numtest,numtest>=1?log10((float)numtest):-1, numill);
      fflush(fd);fflush(fd1);

      if(topresent){
	if(numdeaths>=trigger_dths){
	  startclock++;
	}

	if(startclock==presentday){
	  printf("%d, %d: numinf=%d, numdeaths=%d, \n", r, m, numinf, numdeaths);
	  avinfs+=numinf;avdths+=numdeaths;
	  break;
	}
      }
    }
    fprintf(fd,"\n");fprintf(fd1,"\n");
    for(i=0;i<MAXINFS;i++){//free memory
      if(inflist[i]==1){
	die(infs[i], &numcurinf);//deallocate
      }
    }
  }
  printf("avinfs=%.4f, avdeaths=%.4f\n", avinfs/((double)num_runs), avdths/((double)num_runs));

  free_infar(infs, 0, MAXINFS-1);
  free((char *) P);
  fclose(fd);fclose(fd1);

  return 0;
}

