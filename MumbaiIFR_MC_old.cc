// This is the first version of MumbaiIFR_MC.cc used in the first
// version of the preprint at 
// https://www.medrxiv.org/content/10.1101/2021.04.08.21255101v1

// Generate IFR estimates for Mumbai using Monte Carlo experiments
// Distributions are put on the following variables
// COVID-19 fatalities during 2020
// delay between prevalence & fatality estimate in IFR calculations
// slum prior infection at the time of the 1st survey
// nonslum prior infection at the time of the 1st survey
// fraction of fatalities from the slums prior to the first serosurvey
// fraction of fatalities from the slums after the first serosurvey
// variation in slum naive IFR after the first serosurvey
// variation in nonslum naive IFR after the first serosurvey

// Compile with g++ using "g++ -lm -std=gnu++11 MumbaiIFR_MC_old.cc -o MumbaiIFR_MC_old"


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
  double slumprevmin=0.45, slumprevmax=0.62; // min and max values of "slumprev"
  double nonslumprev; // nonslum prior infection at the time of the 1st survey
  double nonslumprevmin=0.12, nonslumprevmax=0.2; // min and max values of "nonslumprev"
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
  double slumIFRmid, nonslumIFRmid;
  // slum and nonslum naive IFR values post first serosurvey
  double slumIFR, nonslumIFR;
  // Variation allowed in "slumIFR" and "nonslumIFR" relative to "slumIFRmid" and "nonslumIFRmid"
  // (0.1 = 10%)
  double slumIFRvar=0.1, nonslumIFRvar=0.1;
  // (sero)prevalence in the slums, nonslums and city as a whole by the end of 2020
  double slumprevfinal, nonslumprevfinal, totprevfinal;
  // IFR for the city as a whole during 2020
  // (not divided into slum and nonslum IFR as we have no idea how undercounting varies between slums and nonslums)
  double slumIFRfinal, nonslumIFRfinal, IFRfinal;
  // csv file with data from experiments
  FILE *fd=fopen("MumbaiIFRtmp.csv", "w");
  FILE *fd1=fopen("Mumbaiprevhist.csv", "w");//for histograms
  FILE *fd2=fopen("MumbaiIFRhist.csv", "w");//for histograms

  int i,c;

  // data for histograms on seroprevalence and IFR
  int histoprev[20], histoIFR[20];

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
    histoprev[i]=0;histoIFR[i]=0;
  }

  //2020 slum and nonslum populations inferred from 2011 data
  slum2020=slum2011/tot2011*tot2020;
  nonslum2020=nonslum2011/tot2011*tot2020;

  // Header in csv file
  fprintf(fd, "slumprevfinal, nonslumprevfinal, totprevfinal, IFRfinal, slum naive (2020), nonslum naive (2020)\n");

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
    slumIFRmid=P1_slumdeathsfrac*initf[delay]/(slumprev*slum2020 + P1_slumdeathsfrac*initf[delay]);
    nonslumIFRmid=(1-P1_slumdeathsfrac)*initf[delay]/(nonslumprev*nonslum2020 + (1-P1_slumdeathsfrac)*initf[delay]);
    //fprintf(stderr, "slumIFRmid=%.5f, nonslumIFRmid=%.5f\n", slumIFRmid, nonslumIFRmid);

    // choose slum naive IFR post survey
    slumIFR=unif((1.0-slumIFRvar)*slumIFRmid, (1.0+slumIFRvar)*slumIFRmid, generator);

    // choose nonslum naive IFR post survey
    nonslumIFR=unif((1.0-nonslumIFRvar)*nonslumIFRmid, (1.0+nonslumIFRvar)*nonslumIFRmid, generator);

    // total fatalities post survey ("Period 2")
    newdeaths=finalf[delaymin]-initf[delaymin];

    // slum, nonslum, and city-wide infection rate by the end of 2020 (%)
    // final infection rate = (Period 1 infection + added infections)/total population
    slumprevfinal=(slumprev*slum2020 + postwave_slumfrac*newdeaths/slumIFR)/slum2020*100.0;
    nonslumprevfinal=(nonslumprev*nonslum2020 + (1-postwave_slumfrac)*newdeaths/nonslumIFR)/nonslum2020*100.0;
    totprevfinal=(slumprev*slum2020+nonslumprev*nonslum2020+(1.0-postwave_slumfrac)*newdeaths/nonslumIFR+postwave_slumfrac*newdeaths/slumIFR)/tot2020*100.0;

    // IFR at the end of 2020 (%) ["fatalities" in the denominator
    // makes a marginal difference but added for completeness]
    IFRfinal=fatalities/(totprevfinal*tot2020/100.0 + fatalities)*100.0;

    // print
    fprintf(fd, "%.4e, %.4e, %.4e, %.4e, %.4e, %.4e\n", slumprevfinal, nonslumprevfinal, totprevfinal, IFRfinal, (P1_slumdeathsfrac*initf[delay] + postwave_slumfrac*newdeaths)/slumprevfinal*10000.0/slum2020, ((1-P1_slumdeathsfrac)*initf[delay] + (1-postwave_slumfrac)*newdeaths)/nonslumprevfinal*10000.0/nonslum2020);

    //The histograms
    c=int((totprevfinal - 46.0)/2.0);
    if(c>0 && c<=14)
      histoprev[c]++;
    else if(c<=0)
      histoprev[0]++;
    else if(c>14)
      histoprev[15]++;

    c=int((100.0*IFRfinal - 10.0)/2.0);
    if(c>0 && c<=14)
      histoIFR[c]++;
    else if(c<=0)
      histoIFR[0]++;
    else if(c>14)
      histoIFR[15]++;

  }
  fprintf(fd, "\n");

  // print out the histograms
  fprintf(fd1, "<=48,%.1f\n",(double)(histoprev[0])/numexp*100.0);
  for(i=1;i<14;i++){
    fprintf(fd1, "%d-%d,%.1f\n", 2*i+46, 2*i+48,(double)(histoprev[i])/numexp*100.0);
  }
  fprintf(fd1, ">74,%.1f\n", (double)(histoprev[15])/numexp*100.0);

  fprintf(fd1, "\n");

  fprintf(fd2, "<=0.12,%.1f\n",(double)(histoIFR[0])/numexp*100.0);
  for(i=1;i<14;i++){
    fprintf(fd2, "%.2f-%.2f,%.1f\n", (double)(2*i+10)/100.0, (double)(2*i+12)/100.0,(double)(histoIFR[i])/numexp*100.0);
  }
  fprintf(fd2, ">0.38,%.1f\n", (double)(histoIFR[15])/numexp*100.0);


  fclose(fd);fclose(fd1);fclose(fd2);

  return 0;
}
