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
   Contains the function prototypes for the three modes available in main.cpp

   Jonathan Petticrew, University of Sheffield, 2017.
 */


#ifndef MODEL_H
#define MODEL_H

// device_properties() is defined in device_properties.cpp
void device_properties(int material);

// drift_velocity() is defined in drift_velocity.cpp
void drift_velocity(int material);

// ii_coef() is defined in ii_coef.cpp
void ii_coef(int material);

#endif
