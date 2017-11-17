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
histogram_class.cpp contains the function declarations for the Gaussian Histogram Fitter
See  histogram.h for the class definition
1) histogram(double* data, int size);
	Pass the 1d data array and the number of elements it contains.
	The class will ask for a user input for histogram bin width and output the histgram data and fit to Hist.txt
2) histogram(double* data, int size, double binsize);
	Pass the 1d data array, the number of elements the array contains, and the histogram bin width.
	The class will output  the histgram data and fit to Hist.txt
3) histogram(double* data, int size, char* fname);
	Pass the 1d data array, the number of elements the array contains, and a file name to output the data and fit.
	The class will ask for a user input for histogram bin width and output the fit and data to the passed file name.
4) histogram(double* data, int size, double binsize, char* fname);
	Pass the 1d data array, the number of elements the array contains, the histogram bin width, and a file name to output the data and fit.
	The class will utput the fit and data to the passed file name.
	
	Jonathan Petticrew, University of Sheffield, 2017.
*/

#include "histogram.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

//constructor for histogram class Ask for Bin Size
histogram::histogram(double* data, int size):size(size){
	dataset = new double[size];
	int i;
	for(i=0;i<size;i++){
		dataset[i]=data[i];
	}
	mean=0;
	meansquare=0;
	max_min();
	printf("Enter Bin Width:");
	scanf("%lf",&binsize);
	binner();
	fit();
	print();
};

//constructor for histogram class passed Bin Size
histogram::histogram(double* data, int size, double binsize):size(size), binsize(binsize){
	dataset = new double[size];
	int i;
	for(i=0;i<size;i++){
		dataset[i]=data[i];
	}
	mean=0;
	meansquare=0;
	max_min();
	binner();
	fit();
	print();
};

//constructor for histogram class Ask for Bin Size and passes string for printing Hist to file.
histogram::histogram(double* data, int size, char* fname):size(size){
	filename=new char[strlen(fname)+1];
	strcpy(filename,fname);
	dataset = new double[size];
	int i;
	for(i=0;i<size;i++){
		dataset[i]=data[i];
	}
	mean=0;
	meansquare=0;
	max_min();
	printf("Enter Bin Width:");
	scanf("%lf",&binsize);
	binner();
	fit();
	print_custom();
};

//constructor for histogram class passed Bin Size and passes string for printing Hist to file.
histogram::histogram(double* data, int size, double binsize, char* fname):size(size), binsize(binsize){
	filename=new char[strlen(fname)+1];
	strcpy(filename,fname);
	dataset = new double[size];
	int i;
	for(i=0;i<size;i++){
		dataset[i]=data[i];
	}
	mean=0;
	meansquare=0;
	max_min();
	binner();
	fit();
	print_custom();	
};

//destructor for histogram class
histogram::~histogram(){
	delete[] dataset;
	delete[] binEdges;
	delete[] binCenters;
	delete[] binValues;
	if(charstored!=0) delete[] filename;
};

//finds maximum and minimum in dataset
void histogram::max_min(){
	max=dataset[0];
	min=dataset[0];
	int i;
	for(i=1;i<size;i++){
		if(dataset[i]>max) max=dataset[i];
		else if(dataset[i]<min) min=dataset[i];
	}
};

//Bins the DataSet and calculates mean and standard deviation
void histogram::binner(){
	bins=(max-min)/binsize+3;
	binEdges = new double[bins+1];
	binCenters = new double[bins];
	binValues = new double[bins];
	int i, j;
	for(i=0;i<(bins+1);i++){
		binEdges[i]=(min-binsize)+i*binsize;
	}
	for(i=0;i<bins;i++){
		binCenters[i]=(binEdges[i+1]-binEdges[i])/2 +binEdges[i];
		binValues[i]=0;
	}
	
	for(j=0;j<size;j++){
		for(i=0;i<bins;i++){
			if(dataset[j]>=binEdges[i] && dataset[j]<binEdges[i+1]){
				binValues[i]++;
				mean+=binCenters[i]/size;
			}
		}
	}
	for(i=0;i<bins;i++){
		meansquare+=binValues[i]*(binCenters[i]-mean)*(binCenters[i]-mean);
	}
	sdev=sqrt((meansquare-mean)/(size-1));
};

//Fits Gaussian Probability Density Function in the form f(x) = a exp(-((x-b)/c)^2)
void histogram::fit(){
	a=1/(sdev*sqrt(2*3.141592653589793));
	b=mean;
	c=sqrt(2)*sdev;
	FWHM = 2.3548*sdev;
};

//Prints Histogram Data to Hist.txt
void histogram::print(){
	charstored=0;
	FILE *binout;
	binout=fopen("Hist.txt","w");
	fprintf(binout,"f(x)= %lf * exp(-((x-%lf)/%lf)^2)\n",a,b,c);
	int i;
	for(i=0;i<bins;i++){
		fprintf(binout,"%lf, %lf\n",binCenters[i], binValues[i]);
	}
	fclose(binout);
};

//Prints Histogram Data to custom file
void histogram::print_custom(){
	charstored=1; //tells destructor to delete filename
	FILE *binout;
	binout=fopen(filename,"w");
	fprintf(binout,"f(x)= %lf * exp(-((x-%lf)/%lf)^2)\n",a,b,c);
	int i;
	for(i=0;i<bins;i++){
		fprintf(binout,"%lf, %lf\n",binCenters[i], binValues[i]);
	}
	fclose(binout);
};

//Displays Histogram fit in terminal
void histogram::show_fit(){
	printf("f(x)= %lf * exp(-((x-%lf)/%lf)^2)\n",a,b,c);
};
