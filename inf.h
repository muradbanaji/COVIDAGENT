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

#include <stdio.h>
#include <string.h> // for strcmp


class inf{

 protected:

 public:
  int age;  
  int num;
  int numtoinf;
  int ill;
  int quar;
  int inftimes[31];
  int infnums[21];// number to infect at each time

  //constructors, etc

  inf(int orgnum, int P[], int maxP);
  void setinftimes(int rmin, int rmax);

};
