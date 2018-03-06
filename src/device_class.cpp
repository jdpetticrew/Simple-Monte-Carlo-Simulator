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
device_class.cpp contains the class definition of the device class for the SMC.

The device class contains a Poisson solver and a linear interpolater to claculate the
electric field properties of semiconductor devices.

device.cpp contains the class implimentation

Jonathan Petticrew, University of Sheffield, 2017.
*/

#include "device.h"
#include <math.h>
#include <stdio.h>
#include <conio.h>
//Constructor, sets a few variables and calls read() to read in and populate the doping profile.
//PUBLIC
device::device(SMC *input):constants(input){
    q=constants->Get_q();
    Vbi=constants->Get_Vbi();
	  int d1=constants->Get_die();
	  die=d1*8.85e-12;
	  read();
};

//Efield_at_x returns the electric field to main for a given xposition inside
//the electric field profile. If outside the electric field profile it will return 0.
//PUBLIC
double device::Efield_at_x(double xpos){
    double Field;

    int i,j;
    for(i=i_min;i<i_max;i++)
    {
		if (xpos>=efield_x[i])
		{
			j=i+1;
			if(xpos<efield_x[j])
			{
						Field = LinearInterpolate(efield_e[i], efield_e[j], efield_x[i],efield_x[j],xpos);
						return Field;
			}
		}
	}
	return 0;
};
//Get function to return the depletion width.
//PUBLIC
double device::Get_width(){return efield_x[i_max]-efield_x[i_min];
};
double device::Get_xmin(){return efield_x[i_min];
};
double device::Get_xmax(){return efield_x[i_max];
};
//Linear Interpolator used by Efield_at_x
//PRIVATE
double device::LinearInterpolate(double y1, double y2, double x1, double x2, double x){
	double mu, ans;
	mu = (x-x1)/(x2-x1);
	ans= (y1*(1-mu)+y2*mu);
	return ans;
};


//read, populates the doping profile, declares the size of arrays that depend on the doping profile
//Shifts the doping profile read in as cm-3, um to m-3, m and populates the check voltages.
//The check voltages are the points that the device fully depetes regions, i.e. 0 is 2 regions, 1 is 3 regions etc.
//PRIVATE
void device::read(){
	NumLayers=0;
    FILE* doping;
	int inputkey;
	if ((doping=fopen("doping_profile.txt","r"))==NULL){
		printf("Error: Can't open doping profile\n");
		printf("Press space to exit\n");
        while((inputkey=_getch())==0);
	}
	double f,g;
	while(fscanf(doping,"%lf,%lf\n",&f,&g)>0){
		++NumLayers;
	}
	fclose(doping);
	efield_x= new double[NumLayers+1];
	efield_e=new double[NumLayers+1];
	int i;
	for(i=0;i<NumLayers+1;i++){
		efield_x[i]=0;
		efield_e[i]=0;
	}
	N=new double[NumLayers];
	w=new double[NumLayers];
	doping=fopen("doping_profile.txt","r");
	int z;
	for(z=0;z<NumLayers;z++){
		fscanf(doping,"%lf,%lf\n",&N[z],&w[z]);
	}
	fclose(doping);
	for(z=0;z<NumLayers;z++){
	    N[z]=N[z]*1e6;
		  w[z]=w[z]*1e-6;
	}
};
//profiler is the N_Layer electric field solver. Calculates field for voltagein as a
//N-I-P but is required as P-I-N so is then flipped and shifted so the first xposition is 0.
//PUBLIC
void device::profiler(double voltagein){
	double voltage;
	voltage=voltagein+Vbi;
	int i,j;
  int pn=depletionlookup()+1;
  double Vsum=0;
  int err=0;
  double xold=0;
  double xoldest=0;
  int endpoint,depleted;
  while (Vsum<voltage) {
      if (err==0) {
        if (pn>0) {
          pn=pn-1;
        }
        else {
          printf("ERROR, stepped back outside 1st region\n");
        }
      }
      for(i=0;i<(NumLayers+1);i++){
        efield_x[i]=0;
        efield_e[i]=0;
      }
      if (err==0) {
        efield_x[pn+1]=w[pn];
      } else {
        efield_x[pn+1]=xold;
      }
      efield_e[pn+1]=efield_x[pn+1]*q*N[pn]/die;
      j=pn+1;
      depleted=0;
      err=0;
      double ereg;
      while (depleted==0 && err==0) {
        ereg=efield_e[j]+q*N[j]*w[j]/die;
        if (ereg*efield_e[j]>0) {
          efield_e[j+1]=ereg;
          efield_x[j+1]=efield_x[j]+w[j];
          j++;
          if (j>NumLayers-1) {
            err=1;
            xold=efield_x[pn+1]/2;
          }
        }
        else{
          efield_e[j+1]=0;
          efield_x[j+1]=efield_x[j]+fabs((0-efield_e[j])/((efield_e[j]-ereg)/w[j]));
          depleted=1;
        }
      }

      if (err==0) {
        Vsum=0;
        for(i=0;i<NumLayers;i++){
          Vsum=Vsum+0.5*(efield_x[i+1]-efield_x[i])*(efield_e[i]+efield_e[i+1]);
          Vsum=fabs(Vsum);
          //printf("%g\n",Vsum);
        }
      }
      else{
        Vsum=0;
      }
  }
  int solve=0;
  int firstloop=0;
  double xtest=efield_x[pn+1];
  while (solve==0) {
    if (Vsum>voltage) {
      if (firstloop==0) {
        xold=xtest;
        xtest=xtest/2;
        firstloop=1;
      } else {
        xoldest=xold;
        xold=xtest;
        xtest=xold-fabs(xoldest-xold)/2;
      }
    } else {
      if (firstloop==0) {
        xold=xtest;
        xtest=3*xtest/2;
        firstloop=1;
      } else {
        xoldest=xold;
        xold=xtest;
        xtest=xold+fabs(xoldest-xold)/2;
      }
    }
    for(i=0;i<(NumLayers+1);i++){
      efield_x[i]=0;
      efield_e[i]=0;
    }
    efield_x[pn+1]=xtest;
    efield_e[pn+1]=efield_x[pn+1]*q*N[pn]/die;
    int j=pn+1;
    depleted=0;
    while (depleted==0) {
      double ereg=efield_e[j]+q*N[j]*w[j]/die;
      if (ereg*efield_e[j]>0) {
        efield_e[j+1]=ereg;
        efield_x[j+1]=efield_x[j]+w[j];
        j++;
        if (j>NumLayers) {
            printf("Error in field solver\n");
        }
      }
      else{
        efield_e[j+1]=0;
        efield_x[j+1]=efield_x[j]+fabs((0-efield_e[j])/((efield_e[j]-ereg)/w[j]));
        depleted=1;
        endpoint=j+1;
      }
    }
    Vsum=0;
    for(i=0;i<NumLayers;i++){
      Vsum=Vsum+0.5*(efield_x[i+1]-efield_x[i])*(efield_e[i]+efield_e[i+1]);
      Vsum=fabs(Vsum);
      //printf("%g\n",Vsum);
    }
    if (fabs(Vsum-voltage)<0.0001) {
      solve=1;
      //printf("I Solve\n");
    }
  }

  efield_x[pn]=0-efield_x[pn+1];
  for(i=pn+1;i<NumLayers+1;i++){
    efield_x[i]=efield_x[i]+efield_x[pn];
  }
  for(i=0;i<pn+1;i++){
    for(j=0;j<NumLayers+1;j++){
      efield_x[j]=efield_x[j]+w[i];
    }
  }
  for(i=pn-1;i>-1;i--){
    efield_x[i]=efield_x[pn];
  }
  for(i=endpoint+1;i<NumLayers+1;i++){
    efield_x[i]=efield_x[endpoint];
  }
  i_min=pn;
  i_max=endpoint;
  FILE *efield;
  /*efield=fopen("e_out.txt","w");
  for(i=0;i<NumLayers+1;i++){
    fprintf(efield, "%g %g\n",efield_x[i],efield_e[i]);
  }*/

};
//depletionlookup finds the pn junction
int device::depletionlookup(){
  int i;
  for(i=0;i<(NumLayers-1);i++){
    int x;
    x=N[i]*N[i+1];
    if(x<=0){
      FILE *h;
      h=fopen("dep.txt","w");
      fprintf(h,"%d\n",i);
      return i;
    }
  }
};
