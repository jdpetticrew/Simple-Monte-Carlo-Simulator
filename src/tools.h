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
tools.h contains the class definition of the tools class for the SMC
The tools class calculates the scattering rates and probabilities for the SMC

tools_class.cpp contains the class implimentation

Jonathan Petticrew, University of Sheffield, 2017.
*/
#ifndef TOOLS_H
#define TOOLS_H
#include "SMC.h"


const int PB_DIM2_SZ=20000;
const int PB_DIM1_SZ=3;
class tools{
      private:
          SMC *constants;
          double pb[PB_DIM1_SZ][PB_DIM2_SZ];//array of probabilities
          double pb2[PB_DIM1_SZ][PB_DIM2_SZ];//array of probabilities
          void e_rate();
          void h_rate();
          double rtotal;
          double rtotal2;
          double my_pow(double base, double exponent);
      public:
          tools(SMC *input);
          int scattering_probability();
          double Get_rtotal();
          double Get_rtotal2();
          double Get_pb(int i, int j);
          double Get_pb2(int i, int j);
};
#endif
