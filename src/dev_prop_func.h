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
dev_prop_func.h contains the function prototypes for functions that are uniqe to the device_properties mode
in device_properties.cpp

functions are implimented in dev_prop_func.cpp

Jonathan Petticrew, University of Sheffield, 2017.
*/

#ifndef DEV_PROP_FUNC_H
#define DEV_PROP_FUNC_H

//Counts the number of bias in the bias_input.txt
int biascounter();

//Reads in User input for time slices
int timesliceread();

//Reads in User input for injection condition
int usDeviceread();

//Reads in user input for simulation time limit
double simulationtimeread();

//Reads in user input for number of trials per voltage
int trialsread();

//Does some post processing at the end of device_properties()
void postprocess(double Vsim[],double simtime, int voltages);

#endif
