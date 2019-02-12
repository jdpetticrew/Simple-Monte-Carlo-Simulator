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
   SMC_class.cpp contains the class implimentation for the SMC Class
   The SMC class contains all the SMC parameters.

   SMC.h contains the definition
   Functions with Get_ are all getter functions.

   Jonathan Petticrew, University of Sheffield, 2017.
 */
#include "SMC.h"
#include <math.h>
#include <stdio.h>
//Constructor sets non-material specific parameters
SMC::SMC(){
	q=1.6e-19;
	hbar=1.0545716818e-34;
	K=1.380622e-23; //boltzmann's constant
	T=300;
	free_mass=9.1093821545e-31;

};

double SMC::Get_N(){
	return N;
};

double SMC::Get_e_meanpath(){
	return e_meanpath;
};

double SMC::Get_h_meanpath(){
	return h_meanpath;
};

double SMC::Get_q(){
	return q;
};

double SMC::Get_e_Eth(){
	return e_Eth;
};

double SMC::Get_h_Eth(){
	return h_Eth;
};

double SMC::Get_hw(){
	return hw;
};

double SMC::Get_e_mass(){
	return e_mass;
};

double SMC::Get_h_mass(){
	return h_mass;
};

double SMC::Get_e_Cii(){
	return e_Cii;
};

double SMC::Get_h_Cii(){
	return h_Cii;
};

double SMC::Get_e_gamma(){
	return e_gamma;
};

double SMC::Get_h_gamma(){
	return h_gamma;
};

int SMC::Get_NUMPOINTS(){
	return NUMPOINTS;
};

double SMC::Get_hbar(){
	return hbar;
};

double SMC::Get_Emax(){
	return Emax;
};

double SMC::Get_Vbi(){
	return Vbi;
};

double SMC::Get_die(){
	return die;
};

//Sets the material specific parameters
void SMC::mat(int x){

	if (x == 1) { // Silicon Parameters
		e_mass=(0.6*free_mass);
		h_mass=(0.9*free_mass);
		e_meanpath=98e-10;
		h_meanpath=68e-10;
		e_Eth=(1.2*q);
		h_Eth=(1.5*q);
		e_Cii=2e12;
		h_Cii=4.4e12;
		e_gamma=3.5;
		h_gamma=3.5;
		hw=(0.063*q);
		MAX_eV=6;
		Vbi=1;
		die=11.90;
	}
	else if (x == 2) {//GaAd
		e_mass=(0.5*free_mass);
		h_mass=(0.5*free_mass);
		e_meanpath=50.4e-10;
		h_meanpath=47.6e-10;
		e_Eth=(1.75*q);
		h_Eth=(1.75*q);
		e_Cii=40e12;
		h_Cii=30e12;
		e_gamma=4;
		h_gamma=4;
		hw=(0.029*q);
		MAX_eV=8.75;
		Vbi=1.2;
		die=12.9;
	}
	else if (x == 3) {//InGaP
		e_mass=(0.7*free_mass);
		h_mass=(0.7*free_mass);
		e_meanpath=55.7e-10;
		h_meanpath=58.2e-10;
		e_Eth=(2.11*q);
		h_Eth=(2.11*q);
		e_Cii=8e12;
		h_Cii=8e12;
		e_gamma=2.3;
		h_gamma=2.3;
		hw=(0.037*q);
		MAX_eV=6;
		Vbi=1.8;
		die=11.8;
	}

	Emax=MAX_eV*q;
	NUMPOINTS = (int)(MAX_eV*1000);
	N=(1/(exp(hw/(K*T))-1));
};
