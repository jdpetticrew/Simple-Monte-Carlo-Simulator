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
   SMC.h contains the class definition for the SMC Class
   The SMC class contains all the SMC parameters.

   SMC_class.cpp contains the implimentation

   Jonathan Petticrew, University of Sheffield, 2017.
 */

#ifndef SMC_H
#define SMC_H
class SMC {
private:
	double q;
	double MAX_eV;   // Maximum energy in electron volts
	double Emax;   // Maximum energy in joules
	int NUMPOINTS;   //number of points
	double hbar;
	double K;    //boltzmann's constant
	double T;   //temperature
	double hw;
	double N;
	double free_mass;
	double e_mass;
	double h_mass;
	double e_meanpath;
	double h_meanpath;
	double e_Eth;
	double h_Eth;
	double e_Cii;
	double h_Cii;
	double e_gamma;
	double h_gamma;
	double Vbi;
	double die;
public:
	SMC();
	double Get_N();
	double Get_e_meanpath();
	double Get_q();
	double Get_e_Eth();
	double Get_h_Eth();
	double Get_hw();
	double Get_e_mass();
	double Get_h_mass();
	double Get_e_Cii();
	double Get_h_Cii();
	double Get_e_gamma();
	double Get_h_gamma();
	double Get_h_meanpath();
	int Get_NUMPOINTS();
	double Get_hbar();
	double Get_Emax();
	double Get_permitivity();
	double Get_Vbi();
	void mat(int x);
	double Get_die();
};

#endif
