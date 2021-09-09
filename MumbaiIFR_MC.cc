// This is the latest version of this program associated with the 
// revised preprint at https://www.medrxiv.org/content/10.1101/2021.04.08.21255101
// The program and preprint were updated to make use of monthly 
// all-cause mortality data. In particular, a new aim is to 
// separately estimate IFR and the mortality impact in slum 
// and nonslum areas of the city

// Rather than estimating what fraction of excess deaths were 
// from COVID-19, we simply compute estimates based on official 
// fatalities and all-cause excess deaths separately.

// Distributions are put on the following variables
// COVID-19 fatalities during 2020 [no longer used but kept for legacy]
// delay between prevalence & fatality estimate in IFR calculations
// slum prior infection at the time of the 1st survey
// nonslum prior infection at the time of the 1st survey
// fraction of fatalities from the slums prior to the first serosurvey
// fraction of fatalities from the slums after the first serosurvey
// variation in slum naive IFR after the first serosurvey
// variation in nonslum naive IFR after the first serosurvey
// reduction in excess-deaths based IFR during period 2 in slums and nonslums
// baseline mortality during Periods 1 and 2
// the change in death registration coverage during 2020

// Compile with g++ using "g++ -lm -std=gnu++11 MumbaiIFR_MC.cc -o MumbaiIFR_MC"


#include <iostream>
#include <random>

// Real uniform distribution on [lend,rend] 
double unif(double lend, double rend, std::default_random_engine & generator)
{
  static std::uniform_real_distribution<double> dist1;
  return dist1(generator, std::uniform_real_distribution<double>::param_type(lend, rend));
}

// Integer uniform distribution on [lend,rend] 
int unifi(int lend, int rend, std::default_random_engine & generator)
{
  static std::uniform_int_distribution<int> dist1;
  return dist1(generator, std::uniform_int_distribution<int>::param_type(lend, rend));
}

int main()
{
  const int numexp=100000; // number of experiments
  int totsteps=0;
  // delay between prevalence & fatality estimate in IFR calculations
  int delay; 
  int delaymin=0, delaymax=27; // min and max values of "delay"
  double slumprev; // slum prior infection at the time of the 1st survey
  double slumprevmin=0.50, slumprevmax=0.62; // min and max values of "slumprev"
  double nonslumprev; // nonslum prior infection at the time of the 1st survey
  double nonslumprevmin=0.13, nonslumprevmax=0.2; // min and max values of "nonslumprev"
  int fatalities; // Estimated COVID-19 fatalities during 2020
  int fatalitiesmin=11000, fatalitiesmax=25000; // min and max values of "fatalities"
  // 2011 slum, nonslum and total populations in the city
  // from https://portal.mcgm.gov.in/irj/go/km/docs/documents/MCGM%20Department%20List/Public%20Health%20Department/Docs/Census%20FAQ%20%26%20Answer.pdf
  double slum2011=6534460, nonslum2011=5907913, tot2011=12442373;
  // 2020 slum, nonslum and total populations in the city
  // total from https://portal.mcgm.gov.in/irj/portal/anonymous/qlvitalstatsreport?guest\_user=english
  // slum and nonslum fractions assumed to be as in 2011
  double slum2020, nonslum2020, tot2020=12875213;
  // The fraction of reported COVID-19 fatalities after the second serosurvey but during 2020 which occurred in the slums
  double postwave_slumfrac;
  // min and max values of "postwave_slumfrac"
  double postwave_slumfracmin=0.05, postwave_slumfracmax=0.35;
  // fraction of recorded fatalities in the slums during Period 1, and variation in this (0.05=5%)
  // Assumes 52\% of recorded fatalities came from the slums by the time of the serosurvey
  // The 52% value is based on prevalence & naive IFR estimates in Malani et al
  // A. Malani, D. Shah, G. Kang et al. Seroprevalence of SARS-CoV-2 in slums versus non-slums 
  // in Mumbai, India. The Lancet Global Health, 9(2):e110–111, 2021
  // i.e., (0.541*slum2020*0.00076)/(0.541*slum2020*0.00076+0.016*nonslum2020*0.00263) ~ 0.52
  double P1_slumdeathsfrac, P1_slumdeathsfrac_mid=0.52, P1_slumdeathsfrac_var=0.05;
  // total post-serosurvey deaths in 2020 (with delay taken into account)
  int newdeaths;
  // slum and nonslum naive IFR values inferred from 1st serosurvey data
  double wave1slumnIFR, wave1nonslumnIFR;
  // slum and nonslum naive IFR values post first serosurvey
  double wave2slumnIFR, wave2nonslumnIFR;
  // Variation allowed in "wave2slumnIFR" and "wave2nonslumnIFR" relative to "wave1slumnIFR" and "wave1nonslumnIFR"
  // (0.1 = 10%)
  double slumIFRvardown=0.1, nonslumIFRvardown=0.1;
  double slumIFRvarup=0.1, nonslumIFRvarup=0.1;
  // (sero)prevalence in the slums, nonslums and city as a whole by the end of 2020
  double slumprevfinal, nonslumprevfinal, totprevfinal;
  // IFR for the city as a whole during 2020
  // (not divided into slum and nonslum IFR as we have no idea how undercounting varies between slums and nonslums)
  double slumIFRfinal, nonslumIFRfinal, IFRfinal;

  //drop in excess deaths-based IFR in slums and nonslums after the first wave
  //To allow for reduction in overload on hospital system
  double slumeIFRredmin=0.6, slumeIFRredmax=1, slumeIFRred;
  double nonslumeIFRredmin=0.8, nonslumeIFRredmax=1, nonslumeIFRred;

  //drop in registration during pandemic (from an assumed 100%)
  double regdropmin=0.0, regdropmax=0.05, regdrop;

  //baseline registrations (Apr-mid-Jul, mid-Jul-Dec)
  double baseline1mid=24659, baseline2mid=41947;

  //baseline can vary by up to this percentage each side of measured
  double baselinevar=2;
  double baseline1, baseline2;

  //2020 registrations (Apr-mid-Jul, mid-Jul-Dec)
  double pand1=39159, pand2=50200;

  //excess1 = excess at the time of the first serosurvey
  double excess1;
  //excess2 = added excess after first serosurvey
  double excess2;
  
  //infections during first and second periods (n1n2, s1s2 are ratios)
  double sluminfs1, sluminfs2, nonsluminfs1, nonsluminfs2, n1n2, s1s2;

  //total naive slum and nonslum infection fatality rates
  double slumnIFR, nonslumnIFR, nIFR;

  //total excess deaths based slum and nonslum infection fatality rates (wave 1, wave 2 and whole period)
  double wave1slumeIFR, wave1nonslumeIFR;
  double wave2slumeIFR, wave2nonslumeIFR;
  double slumeIFR, nonslumeIFR;
  //In the city as a whole
  double wave1eIFR, eIFR;

  //undercount ratio (ratio of excess deaths based to naive IFR)
  double slumundercount, nonslumundercount, undercount;

  // csv file with data from experiments
  FILE *fd=fopen("MumbaiIFRtmp.csv", "w");
  FILE *fd1=fopen("Mumbaiprevhist.csv", "w");//for histograms
  FILE *fd2=fopen("MumbaiIFRhist.csv", "w");//for histograms

  int i,c;

  int accept=0;

  // data for histograms on seroprevalence and IFR
  int histoprev[20], histoIFR[20], histoeIFR[20];

  std::default_random_engine generator;

  // For random seeding
  int timeint;
  time_t timepoint;

  //recorded fatalities 0 to 27 days after July 8 2020, according to MCGM bulletins
  int initf[28]={5061,5129,5202,5241,5285,5332,5402,5464,5520,5582,5647,5711,5752,5814,5872,5927,5981,6033,6090,6129,6184,6244,6297,6350,6395,6444,6490,6546};
  //recorded fatalities 0 to 27 days after Jan 1 2020, according to MCGM bulletins
  int finalf[28]={11125,11132,11135,11138,11147,11155,11162,11171,11180,11186,11195,11202,11210,11219,11227,11235,11242,11249,11257,11266,11276,11285,11293,11300,11307,11313,11319,11326};

  // Random seeding
  timeint = time(&timepoint); 
  srand(timeint);
  generator.seed(timeint);

  // initialise the histograms
  for(i=0;i<20;i++){
    histoprev[i]=0;histoIFR[i]=0;histoeIFR[i]=0;
  }

  //2020 slum and nonslum populations inferred from 2011 data
  slum2020=slum2011/tot2011*tot2020;
  nonslum2020=nonslum2011/tot2011*tot2020;

  // Header in csv file
  fprintf(fd, "excess deaths, drop in registration (%%), drop in slum eIFR after wave 1 (%%), drop in nonslum eIFR after wave 1 (%%), city-wide infection rate (%%), slum infection rate (%%), nonslum infection rate (%%), city-wide eIFR (%%), slum eIFR (%%), nonslum eIFR (%%), ratio of nonslum to slum eIFR (%%), city-wide wave 1 eIFR (%%), slum wave 1 eIFR (%%), nonslum wave 1 eIFR (%%), city-wide naive IFR (%%), slum naive IFR (%%), nonslum naive IFR (%%), city-wide undercount factor, slum undercount factor, nonslum undercount factor, fraction of excess deaths from the slums (%%), IFRfinal (%%)\n");

  // Run "numexps" experiments
  for(totsteps=0;totsteps<numexp;totsteps++){
    //choose slum seroprevalence from uniform distribution
    slumprev=unif(slumprevmin, slumprevmax, generator);
    //choose nonslum seroprevalence from uniform distribution
    nonslumprev=unif(nonslumprevmin, nonslumprevmax, generator);
    //choose delay from uniform distribution
    delay=unifi(delaymin, delaymax, generator);
    //choose fatalities from uniform distribution
    fatalities=unifi(fatalitiesmin, fatalitiesmax, generator);
    //choose fraction of postwave slum deaths from uniform distribution
    postwave_slumfrac=unif(postwave_slumfracmin, postwave_slumfracmax, generator);
    //choose fraction of recorded Period 1 deaths from the slums
    P1_slumdeathsfrac=unif((1.0-P1_slumdeathsfrac_var)*P1_slumdeathsfrac_mid, (1.0+P1_slumdeathsfrac_var)*P1_slumdeathsfrac_mid, generator);


    // Estimated naive IFR in slums and nonslums at the time of the survey
    // Depends on estimated seroprevalences and the delay
    wave1slumnIFR=P1_slumdeathsfrac*initf[delay]/(slumprev*slum2020 + P1_slumdeathsfrac*initf[delay]);
    wave1nonslumnIFR=(1-P1_slumdeathsfrac)*initf[delay]/(nonslumprev*nonslum2020 + (1-P1_slumdeathsfrac)*initf[delay]);
    //fprintf(stderr, "wave1slumnIFR=%.5f, wave1nonslumnIFR=%.5f\n", wave1slumnIFR, wave1nonslumnIFR);

    // choose slum naive IFR post survey
    wave2slumnIFR=unif((1.0-slumIFRvardown)*wave1slumnIFR, (1.0+slumIFRvarup)*wave1slumnIFR, generator);

    // choose nonslum naive IFR post survey
    wave2nonslumnIFR=unif((1.0-nonslumIFRvardown)*wave1nonslumnIFR, (1.0+nonslumIFRvarup)*wave1nonslumnIFR, generator);

    // total fatalities post survey ("Period 2")
    // corrected this from "delaymin" to "delay"
    newdeaths=finalf[delay]-initf[delay];

    // slum, nonslum, and city-wide infection rate by the end of 2020 (%)
    // final infection rate = (Period 1 infection + added infections)/total population
    slumprevfinal=(slumprev*slum2020 + postwave_slumfrac*newdeaths/wave2slumnIFR)/slum2020*100.0;
    nonslumprevfinal=(nonslumprev*nonslum2020 + (1-postwave_slumfrac)*newdeaths/wave2nonslumnIFR)/nonslum2020*100.0;
    totprevfinal=(slumprev*slum2020+nonslumprev*nonslum2020+(1.0-postwave_slumfrac)*newdeaths/wave2nonslumnIFR+postwave_slumfrac*newdeaths/wave2slumnIFR)/tot2020*100.0;


    //reduction in excess deaths based fatality rates after the first wave period
    slumeIFRred=unif(slumeIFRredmin, slumeIFRredmax, generator);
    nonslumeIFRred=unif(nonslumeIFRredmin, nonslumeIFRredmax, generator);

    //drop in registration baseline
    regdrop=unif(regdropmin, regdropmax, generator);

    //uncertainty in baseline
    baseline1=unif((1-baselinevar/100.0)*baseline1mid,(1+baselinevar/100.0)*baseline1mid,generator);
    baseline2=unif((1-baselinevar/100.0)*baseline2mid,(1+baselinevar/100.0)*baseline2mid,generator);

    //first and second period excess
    excess1=(pand1/(1.0-regdrop)-baseline1);
    excess2=(pand2/(1.0-regdrop)-baseline2);

    //infections at the end of the first and second waves
    sluminfs1=slumprev*slum2020;
    sluminfs2=postwave_slumfrac*newdeaths/wave2slumnIFR;
    nonsluminfs1=nonslumprev*nonslum2020;
    nonsluminfs2=(1-postwave_slumfrac)*newdeaths/wave2nonslumnIFR;

    //excessvector=([excess1],[excess2])
    //infmatrix=([sluminfs1, sluminfs2],[slumeIFRred*sluminfs2, nonslumeIFRred*nonsluminfs2])
    //IFRvector=([wave1slumeIFR],[wave1nonslumeIFR])
    //inverting the equation: excessvector=infmatrix*IFRvector

    //ratios
    n1n2=nonsluminfs1/(nonslumeIFRred*nonsluminfs2);
    s1s2=sluminfs1/(slumeIFRred*sluminfs2);

    //excess deaths based wave 1 IFR in slums and nonslums
    //Solving the equations
    wave1slumeIFR=excess1/(sluminfs1-n1n2*slumeIFRred*sluminfs2) - excess2/(sluminfs1/n1n2-slumeIFRred*sluminfs2);
    wave1nonslumeIFR=excess2/(nonslumeIFRred*nonsluminfs2-nonsluminfs1/s1s2) - excess1/(nonslumeIFRred*nonsluminfs2*s1s2-nonsluminfs1);

    wave2slumeIFR=slumeIFRred*wave1slumeIFR;
    wave2nonslumeIFR=nonslumeIFRred*wave1nonslumeIFR;

    //naive IFR in the slums and nonslum areas
    slumnIFR=(P1_slumdeathsfrac*initf[delay] + postwave_slumfrac*newdeaths)/slumprevfinal*100.0/slum2020;
    nonslumnIFR=((1-P1_slumdeathsfrac)*initf[delay] + (1-postwave_slumfrac)*newdeaths)/nonslumprevfinal*100.0/nonslum2020;
    nIFR=finalf[delay]/(totprevfinal/100.0*tot2020);


    // IFR at the end of 2020 (%) ["fatalities" in the denominator
    // makes a marginal difference but added for completeness]
    IFRfinal=fatalities/(totprevfinal*tot2020/100.0 + fatalities)*100.0;

    //slum and nonslum eIFR over the whole period
    slumeIFR=(wave1slumeIFR*sluminfs1+slumeIFRred*wave1slumeIFR*sluminfs2)/(sluminfs1+sluminfs2);
    nonslumeIFR=(wave1nonslumeIFR*nonsluminfs1+nonslumeIFRred*wave1nonslumeIFR*nonsluminfs2)/(nonsluminfs1+nonsluminfs2);


    if(wave1nonslumeIFR>=wave1nonslumnIFR && wave1slumeIFR>=wave1slumnIFR && wave2nonslumeIFR>=wave2nonslumnIFR && wave2slumeIFR>=wave2slumnIFR && nonslumeIFR>=nonslumnIFR && slumeIFR>=slumnIFR ){//valid simulations
      accept++;

      //should equal excess1+excess2: wave1slumeIFR*sluminfs1+wave2slumeIFR*sluminfs2+wave1nonslumeIFR*nonsluminfs1+wave2nonslumeIFR*nonsluminfs2;



      //wave 1 eIFR and total eIFR
      wave1eIFR=excess1/(sluminfs1+nonsluminfs1);
      eIFR=(excess1+excess2)/(slumprevfinal/100.0*slum2020+nonslumprevfinal/100.0*nonslum2020);
      // print
      fprintf(fd, "%.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e, %.4e\n", excess1+excess2, regdrop*100.0, (1.0-slumeIFRred)*100.0, (1.0-nonslumeIFRred)*100.0, totprevfinal, slumprevfinal, nonslumprevfinal, eIFR*100.0, slumeIFR*100.0, nonslumeIFR*100.0, nonslumeIFR/slumeIFR, wave1eIFR*100.0, wave1slumeIFR*100.0, wave1nonslumeIFR*100.0, nIFR*100.0, slumnIFR*100.0, nonslumnIFR*100.0, eIFR/nIFR, slumeIFR/slumnIFR, nonslumeIFR/nonslumnIFR, (wave1slumeIFR*sluminfs1+wave2slumeIFR*sluminfs2)/(excess1+excess2)*100.0, IFRfinal);

      //The histograms
      c=int((totprevfinal - 44.0)/2.0);
      if(c>0 && c<=14)
	histoprev[c]++;
      else if(c<=0)
	histoprev[0]++;
      else if(c>14)
	histoprev[15]++;

      c=int((10000.0*eIFR - 20)/2.0);
      if(c>0 && c<=14)
	histoeIFR[c]++;
      else if(c<=0)
	histoeIFR[0]++;
      else if(c>14)
	histoeIFR[15]++;
    }

  }
  fprintf(fd, "\n");

  // print out the histograms
  fprintf(fd1, "<=46,%.1f\n",(double)(histoprev[0])/accept*100.0);
  for(i=1;i<14;i++){
    fprintf(fd1, "%d-%d,%.1f\n", 2*i+44, 2*i+46,(double)(histoprev[i])/accept*100.0);
  }
  fprintf(fd1, ">72,%.1f\n", (double)(histoprev[15])/accept*100.0);

  fprintf(fd1, "\n");

  fprintf(fd2, "<=0.22,%.1f\n",(double)(histoeIFR[0])/accept*100.0);
  for(i=1;i<14;i++){
    fprintf(fd2, "%.2f-%.2f,%.1f\n", (double)(2*i+20)/100.0, (double)(2*i+22)/100.0,(double)(histoeIFR[i])/accept*100.0);
  }
  fprintf(fd2, ">0.48,%.1f\n", (double)(histoeIFR[15])/accept*100.0);

  fprintf(stderr, "%d valid experiments out of %d.\n", accept, numexp);

  fclose(fd);fclose(fd1);fclose(fd2);

  return 0;
}
