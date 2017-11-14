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
#include <assert.h>
#include "functions.h"
#include "math.h"
//All these functions are for Get, Set and Zero.
//Assert is being used to pick up array overflow vilations which were causing problems.
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
	 assert(i<Array && "Position array violated");
     position[i]=input;
};
void carrier::Input_Egy(int i, double input){
	assert(i<Array && "Egy array violated");
     Egy[i]=input;
};
void carrier::Input_kxy(int i, double input){
	assert(i<Array && "kxy array violated");
     kxy[i]=input;
};
void carrier::Input_scattering(int i, int input){
	assert(i<Array && "Scattering array violated");
     scattering[i]=input;
};
void carrier::Input_kz(int i, double input){
	assert(i<Array && "kz array violated");
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
	assert(i<Array && "time array violated");
	time[i]=input;
};
double carrier::Get_time(int i){
	return time[i];
};
void carrier::Input_dt(int i, double input){
	assert(i<Array && "dt array violated");
	dt[i]=input;
};
double carrier::Get_dt(int i){
	return dt[i];
};
void carrier::Input_dx(int i, double input){
	assert(i<Array && "dx array violated");
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
	assert(i<Array && "timearray violated");
	timearray[i]=input;
};
int carrier::Get_timearray(int i){
	return timearray[i];
};

//Calculates the new scattering direction and momentums
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

