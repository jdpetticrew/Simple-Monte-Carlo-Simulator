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
    jcheck=0;
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
    for(i=0;i<(i_max);i++)
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
double device::Get_width(){return width;
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
	G=new double[NumLayers];
	Gw=new double[NumLayers];
	Nwe=new double[NumLayers];
	Nw2e=new double[NumLayers];
	for(i=0;i<NumLayers;i++){
		G[i]=q*N[i]/die;
		G[0]=fabs(G[0]);
		Gw[i]=G[i]*w[i];
		Nwe[i]=N[i]*w[i]/die;
		Nw2e[i]=N[i]*w[i]*w[i]/die;
	}

	Vtcheck=new double[NumLayers-1];
	Vtcheck[0]=0.5*(w[1]+Gw[1]/G[0])*Gw[1];
	for(i=1;i<NumLayers-1;i++){
		Vtcheck[i]=vtr(i);
	}
};
//profiler is the N_Layer electric field solver. Calculates field for voltagein as a
//N-I-P but is required as P-I-N so is then flipped and shifted so the first xposition is 0.
//PUBLIC
void device::profiler(double voltagein){
	double voltage;
	voltage=voltagein+Vbi;
	int i,numlayers;
	for(i=0;i<(NumLayers+1);i++){
		efield_x[i]=0;
		efield_e[i]=0;
	}
};
