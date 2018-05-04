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
   functions.h contains the function prototypes for the common functions used in all three modes
   function declerations in functions.cpp

   Jonathan Petticrew, University of Sheffield, 2017.
 */

#ifndef FUNC_H
#define FUNC_H

/* Random Generator Initilization*/
/* Period parameters */
#define Nr 624
#define M 397
#define MATRIX_A 0x9908b0df     /* constant vector a */
#define UPPER_MASK 0x80000000   /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff   /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


static unsigned long mt[Nr];
static int mti = Nr+1;
double _max(double x, double y);
void sgenrand(unsigned long seed);
double genrand();
#endif
