/* Copyright (C) 2020, Murad Banaji
 *
 * This file is part of COVIDAGENT
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

#define MAXAGE 25
#define MAXDISCPROB 120

class inf{

 protected:

 public:
  int age;//days since infection. Updates at start of iteration.
  int num;//index in list
  int numtoinf;//number who will be infected (without mitigation)
  int ill;
  int quar;//in quarantined state?
  int quardt;//quarantine date
  int inftimes[MAXDISCPROB+1];
  int infnums[MAXAGE];// number to infect at each time
  int dth_time;
  int recov_time;
  int sero_time;//time to seroconversion (currently can't be longer than lifespan)
  int lastop_time;// when can it be destroyed?

  //constructors, etc

  inf(int orgnum, int P[], int maxP);//arbitrary distribution
  inf(int orgnum, double alpha, double beta);//gamma distribution
  void setinftimes(int rmin, int rmax);//uniform distribution
  void setinftimes(double alpha, double beta);//gamma distribution

};
