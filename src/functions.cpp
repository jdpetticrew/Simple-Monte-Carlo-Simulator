/* Copyright 2017 Advanced Detector Centre, Department of Electronic and
   Electrical Engineering, University of Sheffield, UK.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.*/

/*
   functions.cpp contains the function declerations for the common functions used in all three modes
   function prototypes in functions.h

   Jonathan Petticrew, University of Sheffield, 2017.
 */
#include "functions.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <conio.h>

double _max(double x,double y)
{
	if(x>y) return x;
	else return y;
};
// seeds the random number generator
void sgenrand(unsigned long seed)
{
	int i;
	for (i=0; i<Nr; i++)
	{    mt[i] = seed & 0xffff0000;
		 seed = 69069*seed+1;
		 mt[i]=(seed & 0xffff0000)>>16;
		 seed = 69069*seed+1;}
	mti=Nr;
};

//calculates the next random number in the seeded mersenne twister.
double genrand(){ //quite a lot of the parameters in this function are #defined in default_include.h
	double ans;
	unsigned long y;
	static unsigned long mag01[2] = {0x0, MATRIX_A};
	if (mti >= Nr) //generate N words at one time
	{    int kk;
		 if (mti == Nr+1) sgenrand(4357); //sets initial seed
		 for (kk=0; kk<(Nr-M); kk++)
		 {    y=(mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		  mt[kk] = mt[kk+M] ^ (y>>1) ^ mag01[y & 0x1];}
		 for (; kk<Nr-1; kk++)
		 {    y=(mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
		  mt[kk] = mt[kk+(M-Nr)] ^ (y>>1) ^ mag01[y & 0x1];}
		 y=(mt[Nr-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		 mt[Nr-1]=mt[M-1]^(y>>1)^mag01[y &0x1];
		 mti=0;}
	y = mt[mti++];
	y ^= TEMPERING_SHIFT_U (y);
	y ^= TEMPERING_SHIFT_S (y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T (y) & TEMPERING_MASK_C;
	y ^= TEMPERING_SHIFT_L (y);
	ans= ( double   )y * 2.3283064365386963e-10;
	while ((ans<=0)||(ans>=1)) ans=genrand();
	return ans;
};
