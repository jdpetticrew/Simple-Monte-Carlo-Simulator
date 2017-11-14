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
device.h contains the class implimentation of the device class for the SMC.

The device class contains a Poisson solver and a linear interpolater to claculate the
electric field properties of semiconductor devices.

device_class.cpp contains the class definition

Jonathan Petticrew, University of Sheffield, 2017.
*/

#ifndef DEVICE_H
#define DEVICE_H
#include "SMC.h"

class device{
      private:
          SMC *constants; 
		  int NumLayers;                  
          double Vbi;         
          double width;        
          double *efield_x;
          double *efield_e;
          double *N;
          double *w;
          double *G;
          double *Gw;
          double *Nwe;
          double *Nw2e;
          double *Vtcheck;
          int i_max;          
          int jcheck;
          double die;
          double q;
          double LinearInterpolate(double y1, double y2, double x1, double x2, double x);
          double vtr(int b);
	      double Wr(double var[], int a, int b);
		  double bmid(int a, int b);
		  double cmid(int k);
		  void read();
      public:
          device(SMC *con);         
          double Efield_at_x(double xpos);//returns the Efield for a given position
          double Get_width();
          void profiler(double voltage);
};
#endif
