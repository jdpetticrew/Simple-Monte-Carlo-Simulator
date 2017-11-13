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

#ifndef DEV_PROP_FUNC_H
#define DEV_PROP_FUNC_H

int biascounter();
int timesliceread();
int usDeviceread();
double simulationtimeread();
int trialsread();
void postprocess(double Vsim[],double simtime, int voltages);

#endif
