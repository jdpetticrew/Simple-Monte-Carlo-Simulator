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
carrier_class.cpp contains the class definition for the carrier class for the SMC
The carrier class contains all the information about the carriers as they travel through the device

carrer.h contains the class implimentation

Jonathan Petticrew, University of Sheffield, 2017.
*/

#include "carrier.h"
#include "functions.h"
#include "math.h"
//All these functions are for Get, Set and Zero.
carrier::carrier(SMC *input):constants(input){
	int i;
	emass= constants->Get_e_mass();
	hmass= constants->Get_h_mass();
	hbar=constants->Get_hbar();
	for(i=1;i<(Array+1);i++)
	{    position[i]=0;
		 Egy[i]=0;
		 kxy[i]=0;
		 kz[i]=0;
		 scattering[i]=0;	
		 dt[i]=0;
		 dx[i]=0;
		 time[i]=0;
		 timearray[i]=0;
	}
};
void carrier::Input_pos(int i, double input){
     position[i]=input;
};
void carrier::Input_Egy(int i, double input){
     Egy[i]=input;
};
void carrier::Input_kxy(int i, double input){
     kxy[i]=input;
};
void carrier::Input_scattering(int i, int input){
     scattering[i]=input;
};
void carrier::Input_kz(int i, double input){
     kz[i]=input;
};
double carrier::Get_pos(int i){return position[i];
};
double carrier::Get_Egy(int i){return Egy[i];
};
double carrier::Get_kxy(int i){return kxy[i];
};
double carrier::Get_kz(int i){return kz[i];
};
int carrier::Get_scattering(int i){return scattering[i];
};
carrier::~carrier(){
};
void carrier::Input_time(int i, double input){
	time[i]=input;
};
double carrier::Get_time(int i){
	return time[i];
};
void carrier::Input_dt(int i, double input){
	dt[i]=input;
};
double carrier::Get_dt(int i){
	return dt[i];
};
void carrier::Input_dx(int i, double input){
	dx[i]=input;
};
double carrier::Get_dx(int i){
	return dx[i];
};
// Resets all the arrays to 0.
void carrier::reset(){
	int i;
	for(i=1;i<(Array+1);i++)
	{    position[i]=0;
		 Egy[i]=0;
		 kxy[i]=0;
		 kz[i]=0;
		 scattering[i]=0;	
		 dt[i]=0;
		 dx[i]=0;
		 time[i]=0;
		 timearray[i]=0;
	}
};
void carrier::Input_timearray(int i, int input){
	timearray[i]=input;
};
int carrier::Get_timearray(int i){
	return timearray[i];
};

//Calculates the new scattering direction and momenta
void carrier::scatter(int i, int j){
	double cos_theta,kf;
    if (j==0 ) kf=2*(emass)*Egy[i]/((hbar)*(hbar));
    else kf=2*(hmass)*Egy[i]/((hbar)*(hbar));
    if(kf>=0){
    	cos_theta=2*genrand()-1;
    	kz[i]=cos_theta*sqrt(kf);
    	kxy[i]=kf*(1-cos_theta*cos_theta);
    }
};
//Generates a new carrier after an impact ionization event
void carrier::generation(int i, double z_pos, double Energy, double timein, double dtin, int timearrayin){
	position[i]=z_pos;
	Egy[i]=Energy;
	time[i]=timein;
	dt[i]=dtin;
	timearray[i]=timearrayin;
};

