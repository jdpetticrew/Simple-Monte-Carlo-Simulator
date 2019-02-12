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
   carrier.h contains the class implimentation for the carrier class for the SMC
   The carrier class contains all the information about the carriers as they travel through the device

   carrer_class.cpp contains the class definition

   Jonathan Petticrew, University of Sheffield, 2017.
 */

#ifndef CARRIER_H
#define CARRIER_H
#define Array 1000000 // Statically allocated max number of electrons or holes to track.
#include "SMC.h"
class carrier {
private:
	double position[Array];
	double Egy[Array];
	double kxy[Array];
	double kz[Array];
	int scattering[Array];
	double time[Array];
	double dt[Array];
	double dx[Array];
	int timearray[Array];
	SMC *constants;
	double hmass;
	double emass;
	double hbar;
public:
	carrier(SMC *con);
	~carrier();
	void Input_pos(int i, double input);
	void Input_scattering(int i, int input);
	void Input_Egy(int i, double input);
	void Input_kxy(int i, double input);
	void Input_kz(int i, double input);
	double Get_pos(int i);
	double Get_Egy(int i);
	double Get_kxy(int i);
	double Get_kz(int i);
	int Get_scattering(int i);
	void Input_time(int i, double input);
	double Get_time(int i);
	void Input_dt(int i, double input);
	double Get_dt(int i);
	void Input_dx(int i, double input);
	double Get_dx(int i);
	void reset();
	void Input_timearray(int i, int input);
	int Get_timearray(int i);
	void scatter(int i, int j);
	void generation(int i, double z_pos, double Egy, double time, double dt, int timearray);
};
#endif
