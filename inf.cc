#include "inf.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h> // random seeding
#include <ctype.h>

// Maximum number of infected individuals
#define MAXINFS 10000000

//Externally declared (bad practice I know!)

int numinf;//number infected (cumulative)
int numcurinf;//number currently infected
int numquar;//number quarantined (cumulative)
int numtest;//number tested (cumulative)
int numill;//number ill (cumulative)
int inflist[MAXINFS]; //For book-keeping free spaces in list


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

//Poisson distribution
double Pois(double lamb, int k){
  return pow(lamb, (double)(k))*exp(-lamb)/(double)(factorial(k));
}

//Geometric distribution including 0
double Geom(double Ep, int k){
  return pow((1.0-1.0/(1.0+Ep)), (double)k)/(1.0+Ep);
}

//Setting the tail to zero in the distributions below lowers the true expected value. Also given the different weights of the tails, it can make the dynamics for the two distributions differ.

//Poisson distribution cast to integers (not rounded) and with the tail given to zero 
int *discPois(double lamb, int max){
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

//Geometric distribution cast to integers (not rounded) and with the tail given to zero
int *discGeom(double lamb, int max){
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

int propill;
int propdeath;

double percill;//percentage who fall (seriously) ill
double percdeath;//percentage of ill who die

int create(inf *infs[], int P[], int maxP){
  int i=firstfreepos(inflist, MAXINFS);
  int j;
  if(i==-1)//no more space
    exit(0);
  infs[i] = new inf(i, P, maxP);
  numinf++;
  numcurinf++;

  for(j=0;j<21;j++){//number to infect at time j
    infs[i]->infnums[j]=0;
  }
  // 1 in every propill will fall ill?
  if(randpercentage(percill)){
    infs[i]->ill=1;
    numill++;
    //fprintf(stderr, "ill=%d\n", infs[i]->ill);
  }
  //set infection times
  infs[i]->setinftimes(3, 14);
  inflist[i]=1;
  return i;

}

void die(inf *a){
  numcurinf--;
  // if(a->ill==1)//if ill
  //   numill--;
  // if(a->quar==1)//quarantined
  //   numquar--;
  inflist[a->num]=0;
  delete a;
  return;
}


int main(int argc, char *argv[]){
  int timeint;
  time_t timepoint;
  int i, j, m, r, cur, rmax=1;//number of runs
  int maxP;
  double lambda=3.5, lold;
  int *P;
  int totdays=150;
  inf **infs;
  int initinfs=10;
  int numdeaths=0;
  int dth=200;
  float dthrate=0.5;//percentage
  int geometric=0;//geometric or poisson?
  int haslockdown=0;
  int lockdownlen=21;
  int lockdown=0;//at 10th death
  int lockdownday=0;
  float lockdowneff=75.0;// as a percentage
  //physical distancing?
  int haspd=1;
  int numpddth=1;//starts at nth death
  float pdeff=30.0;// efficiency. E.g. 40% - removes 2 in 5 contacts
  int pd=0;
  int pdday=0;
  int startclock=0;
  char fname[200];
  double avinfs, avdths;
  FILE *fd;
  int topresent=0;//only simulate to a fixed day, namely "presentday" days after 10th death
  int presentday=18;// days since 10th death
  double quarp=4.0;//chances of getting quarantined as a percentage
  int testdate=12;//Currently assume all tests occur on a particular day in the infection cycle. Only those tested are quarantined. India: 7? UK: 12?
  int numlockdth=200;//The lockdown begins after the death number numlockdeath. UK about 250, India 10
  long totpop=6600000;
  int herd=1;//herd immunity?
  double herdlevel=0;


  if(argc < 5){
    fprintf(stderr, "ERROR: you must enter the death rate (percentage), geometric distribution? (0 for Poisson or 1 for geometric), lambda (R0 float value), haslockdown? (1 for yes). If there is a lockdown, you should enter its lenghth. You may also enter the lockdown efficiency as a percentage. (75 percent is the default.)\n");
    return -1;
  }

  dthrate=atof(argv[1]);
  //dth = atoi(argv[1]);
  geometric=atoi(argv[2]);
  lambda=atof(argv[3]);
  haslockdown=atoi(argv[4]);
  if(haslockdown){
    if(argc < 6){
      fprintf(stderr, "ERROR: you must enter the lockdown length.\n");
      return -1;
    }
    else{
      lockdownlen=atoi(argv[5]);
      if(argc >= 7){
	lockdowneff=atof(argv[6]);
      }
    }
  }
  if(haslockdown)
    sprintf(fname, "data2/dth%.1f_%s_lamb%.2f_lock%d_%.1f.csv", dthrate, geometric==1?"geom":"pois", lambda, lockdownlen, lockdowneff);
  else
    sprintf(fname, "data2/dth%.1f_%s_lamb%.2f_nolock.csv", dthrate, geometric==1?"geom":"pois", lambda);
  fprintf(stderr, "%s\n", fname);

  fd=fopen(fname, "w");

  //How is the number to infect distributed?
  if(geometric){
    maxP=20;
    P=discGeom(lambda, maxP);
  }
  else{//Poisson
    maxP=10;
    P=discPois(lambda, maxP);
  }

  //product of these two is 1/ death rate
  propill=5;
  propdeath=dth/propill;

  percill=20.0;//percentage of people who fall seriously ill
  percdeath=dthrate*100.0/percill;

  timeint = time(&timepoint); /*convert time to an integer */
  srand(timeint);

  for(i=0;i<=maxP;i++){
    fprintf(stderr, "%d\n", P[i]);
  }


//At each time step go through each infected, augment its age, check its numtoinf, and then go through its inftimes, and spawn new if needed, remove it if needed.
  avinfs=0.0;avdths=0.0;
  for(r=0;r<rmax;r++){//Each model run
    startclock=0;
    numinf=0;//cumulative total infections
    numcurinf=0;//current infections
    numquar=0;
    numtest=0;
    numill=0;
    numdeaths=0;
    lockdown=0;
    lockdownday=0;
    infs=infar(0, MAXINFS-1);
    for(i=0;i<MAXINFS;i++){
      inflist[i]=0;
    }

    for(i=0;i<initinfs;i++){
      cur=create(infs, P, maxP);//create
      fprintf(stderr, "infs[%d] (illstate=%d) will infect %d at times:\n",cur, infs[cur]->ill, infs[cur]->numtoinf);
      for(j=0;j<(infs[cur])->numtoinf;j++){
	fprintf(stderr, "   %d\n", infs[cur]->inftimes[j]);
      }
    }
    
    for(m=0;m<totdays;m++){
      if(herd){
	fprintf(stderr, "herdlevel=%.4f\n", herdlevel);
      }

      if(haslockdown){
	if(numdeaths>=numlockdth && lockdownday<lockdownlen){
	  lockdown=1;
	  fprintf(stderr, "in lockdown. pdeff up to 60pc.\n");
	  if(lockdownday==0){
	    totpop/=18;
	    pdeff=60.0;
	  }
	  // else{
	  //   totpop+=10000;
	  // }
	  lockdownday++;
	}
	else
	  lockdown=0;
      }

      if(haspd){
	if(numdeaths>=numpddth){
	  pd=1;
	  fprintf(stderr, "physical distancing.\n");
	}
	else
	  pd=0;
      }


      for(i=0;i<MAXINFS;i++){
	if(inflist[i]==1){
	  (infs[i]->age)++;//age updates at start...?
	  if(infs[i]->age==testdate && randpercentage(quarp)){
	    infs[i]->quar=1;
	    numquar++;
	    numtest++;
	  }
	  if((infs[i]->age>17 && infs[i]->ill==1) || infs[i]->age>20){
	    if((infs[i]->ill)==1){//if ill
	      // if(randnum(propdeath)<1)// 1 in propdeath chance of dying if ill
	      if(randpercentage(percdeath))
		numdeaths++;
	    }
	    die(infs[i]);// die or recover
	  }
	  else if(infs[i]->quar==0){//if not quarantined, then infect
	    if(!pd|| (pd && randpercentage(100.0-pdeff))){//one way of dealing with lockdown
	      if(herd){herdlevel=100.0*((double)numinf/(double)totpop);}
	      if(!herd || (herd && randpercentage(100.0-herdlevel))){
		for(j=0;j<infs[i]->infnums[infs[i]->age];j++)
		  create(infs, P, maxP);
	      }
	    }
	  }
	}
      }

      fprintf(stderr, "%d: numinf=%d, numcurinf=%d, numdeaths=%d, numtest=%d, numill=%d\n", m, numinf, numcurinf, numdeaths, numtest, numill);
      fprintf(fd,"%d,%d,%.4f,%d,%d,%.4f, %d,%d\n", m, numinf, numinf>=1?log10((float)numinf):-1, numcurinf, numdeaths, numdeaths>=1?log10((float)numdeaths):-1,numtest, numill);
      fflush(fd);
      if(topresent){
	if(numdeaths>=10){
	  startclock++;
	}
	if(startclock==presentday){
	  printf("%d, %d: numinf=%d, numdeaths=%d, \n", r, m, numinf, numdeaths);
	  avinfs+=numinf;avdths+=numdeaths;
	  for(i=0;i<MAXINFS;i++){
	    if(inflist[i]==1){
	      die(infs[i]);//deallocate
	    }
	  }
	  break;
	}
      }
    }
    free_infar(infs, 0, MAXINFS-1);
  }
  printf("avinfs=%.4f, avdeaths=%.4f\n", avinfs/((double)rmax), avdths/((double)rmax));

  free((char *) P);
  fclose(fd);

  return 0;
}

