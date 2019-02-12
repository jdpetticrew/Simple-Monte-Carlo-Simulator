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
   dev_prop_func.cpp contains the function definitions for functions that are uniqe to the device_properties mode
   in device_properties.cpp

   functions are prototypes in dev_prop_func.h
   Uses the class histogram in postprocess()
   Jonathan Petticrew, University of Sheffield, 2017.
 */

#include "dev_prop_func.h"
#include "histogram.h"
#include <stdio.h>
#include <tchar.h>
#include <math.h>

//Counts the Number of Bias in bias_input.txt
int biascounter(){
	double voltage;
	FILE *bias;
	if ((bias=fopen("bias_input.txt","r"))==NULL)
	{   printf("Error: bias_input.txt can't be opened'\n");}
	int bias_count=0;
	while(fscanf(bias,"%lf",&voltage)>0) {
		bias_count++;
	}
	fclose(bias);
	return bias_count;
};

//Reads in User inout
int timesliceread(){
	int timeslice;
	printf("How many divisions per transit time: \n");
	scanf("%d",&timeslice);
	return timeslice;
};
//Reads in User inout
int usDeviceread(){
	int usDevice;
	printf("1)Pure Electron, 2)Pure Hole:\n");
	scanf("%d",&usDevice);
	return usDevice;
};
//Reads in User inout
double simulationtimeread(){
	double simulationtime;
	printf("Simulation Time in ps:\n");
	scanf("%lf",&simulationtime);
	simulationtime=simulationtime*1e-12;
	return simulationtime;
};
//Reads in User inout
int trialsread(){
	int Ntrials;
	printf("Number of trials (Default=10000):\n");
	scanf("%d", &Ntrials);
	return Ntrials;
};

//Calculates Gain, Noise, Mean Time (using 0.1ps bin width)
void postprocess(double Vsim[],double simtime, int voltages){
	int numbins=(int)(simtime/0.1e-12);
	int i;
	double *G = new double[voltages];
	double *F = new double[voltages];
	double *T = new double[voltages];
	FILE *results;
	results=fopen("Result_2.txt","w");
	fprintf(results,"Voltage Gain Noise MeanTime(ps)\n");
	for(i=0; i<voltages; i++) {
		double V=Vsim[i];
		FILE *Mout;
		char nameM[]= "gain_out.txt";
		char voltagetb[8];
		snprintf(voltagetb,sizeof(voltagetb),"%g",V);
		int fileM_len = strlen(voltagetb) + strlen(nameM) + 1;
		char *fileM = new char[fileM_len];
		snprintf(fileM,fileM_len,"%s%s",voltagetb,nameM);
		Mout=fopen(fileM,"r");
		delete[] fileM;
		int event, count, count2;
		double scanned, TGain, Gain2,mgain2;
		TGain=0;
		Gain2=0;
		count=0;
		while((fscanf(Mout,"%d %lf %d\n",&event, &scanned, &count2))>0) {
			//scanned=scanned/1.6e-19;
			TGain+=scanned;
			Gain2+=scanned*scanned;
			count++;
		}
		G[i]=TGain/count;
		mgain2=Gain2/count;
		F[i]=mgain2/(G[i]*G[i]);
		fclose(Mout);
		FILE *Tout;
		char nameT[]= "time_to_breakdown.txt";
		int fileT_len = strlen(voltagetb) + strlen(nameT) + 1;
		char *fileT =new char[fileT_len];
		snprintf(fileT,fileT_len,"%s%s",voltagetb,nameT);
		if((Tout=fopen(fileT,"r"))!=NULL) {
			count=0;
			int dump;
			double dump2;
			while(fscanf(Tout,"%d %lf\n",&dump,&dump2)>0) {
				count++;
			}
			if(count>0) {
				double *data = new double[count];
				data = { 0 };
				rewind(Tout);
				count=0;
				while(fscanf(Tout,"%d %lf\n",&dump,&dump2)>0) {
					data[count]=dump2/1e-12;
					count++;
				}
				fclose(Tout);
				delete[] data;
				char nameH[]="Hist.txt";
				int fileH_len = strlen(voltagetb) + strlen(nameH) + 1;
				char *fileH = new char[fileH_len];
				snprintf(fileH, fileH_len,"%s%s",voltagetb,nameH);
				histogram hist(data,count,0.1,fileH);
				delete[] fileH;
				T[i]=hist.Get_Mean();
				fprintf(results,"%lf %lf %lf %lf\n",Vsim[i],G[i],F[i],T[i]);
			}
			else{
				fclose(Tout);
				fprintf(results,"%lf %lf %lf --\n",Vsim[i],G[i],F[i]);
			}
		}
		else fprintf(results,"%lf %lf %lf --\n",Vsim[i],G[i],F[i]);
		delete[] fileT;
	}
	fclose(results);
	delete[] G,F,T;
};
