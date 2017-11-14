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
device_class.cpp contains the class definition of the device class for the SMC.

The device class contains a Poisson solver and a linear interpolater to claculate the
electric field properties of semiconductor devices.

device.cpp contains the class implimentation

Jonathan Petticrew, University of Sheffield, 2017.
*/

#include "device.h"
#include <math.h>
#include <stdio.h>
#include <conio.h>
//Constructor, sets a few variables and calls read() to read in and populate the doping profile.
//PUBLIC
device::device(SMC *input):constants(input){
    q=constants->Get_q();    
    Vbi=constants->Get_Vbi();    
    jcheck=0;
	int d1=constants->Get_die();
	die=d1*8.85e-12;
	read();
};

//Efield_at_x returns the electric field to main for a given xposition inside
//the electric field profile. If outside the electric field profile it will return 0.
//PUBLIC
double device::Efield_at_x(double xpos){
    double Field;
       
    int i,j;
    for(i=0;i<(i_max);i++)
    {
		if (xpos>=efield_x[i])
		{   
			j=i+1;
			if(xpos<efield_x[j])
			{   
						Field = LinearInterpolate(efield_e[i], efield_e[j], efield_x[i],efield_x[j],xpos);
						return Field;
			}
		}				
	}
	return 0;   
};
//Get function to return the depletion width.
//PUBLIC
double device::Get_width(){return width;
};
//Linear Interpolator used by Efield_at_x
//PRIVATE
double device::LinearInterpolate(double y1, double y2, double x1, double x2, double x){
	double mu, ans;
	mu = (x-x1)/(x2-x1);
	ans= (y1*(1-mu)+y2*mu);
	return ans;
};

//vtr is used as part of the electric field calculation for the N_layer solver
//used to calculate check voltages
//PRIVATE
double device::vtr(int b){
	double ans=0;
	ans=0.5*w[b+1]*Gw[b+1]+0.5*w[b]*(Gw[b+1]+Wr(Gw,b,b+1))+0.5*(Wr(Gw,1,b+1)*Wr(Gw,1,b+1))/G[0];
	if (b>1){
		int n;
		for(n=2;n<b+1;n++){
			ans+=0.5*w[n-1]*(Wr(Gw,n,b+1)+Wr(Gw,n-1,b+1));
		}
	}
	
	return ans;
};
//Wr is used as part of the electric field calculation for the N_layer solver
//sums var between a and b
//PRIVATE
double device::Wr(double var[], int a, int b){
	double ans=0;
	int i;
	for(i=a;i<(b+1);i++){
		ans+=var[i];
	}
	return ans;
};
//bmid is used as part of the electric field calculation for the N_Layer solver
//returns part of the quadratic eqn
//PRIVATE
double device::bmid(int a, int b){
	double ans=0;
	int i;
	for (i=a;i<(b+1);i++){
		ans+=w[i]-(Nwe[i]/(N[0]/die));
	}
	return ans;
};
//cmid is used as part of the electric field calculation for the N_Layer solver
//returns part of the quadratic eqn
//PRIVATE
double device::cmid(int k){
	double ans =0;
	int i;
	if(k>1){
		for(i=2;i<(k+1);i++){
			ans+=2*Nwe[i]*Wr(w,1,i-1);
		}
	}
	return ans;
};
//read, populates the doping profile, declares the size of arrays that depend on the doping profile
//Shifts the doping profile read in as cm-3, um to m-3, m and populates the check voltages.
//The check voltages are the points that the device fully depetes regions, i.e. 0 is 2 regions, 1 is 3 regions etc.
//PRIVATE
void device::read(){
	NumLayers=0;
    FILE* doping;
	int inputkey;
	if ((doping=fopen("doping_profile.txt","r"))==NULL){
		printf("Error: Can't open doping profile\n");
		printf("Press space to exit\n");
        while((inputkey=_getch())==0); 
	}
	double f,g;
	while(fscanf(doping,"%lf,%lf\n",&f,&g)>0){
		++NumLayers;
	}
	fclose(doping);
	efield_x= new double[NumLayers+1];
	efield_e=new double[NumLayers+1];
	int i;
	for(i=0;i<NumLayers+1;i++){
		efield_x[i]=0;
		efield_e[i]=0;
	}
	N=new double[NumLayers];
	w=new double[NumLayers];
	doping=fopen("doping_profile.txt","r");
	int z;
	for(z=0;z<NumLayers;z++){
		fscanf(doping,"%lf,%lf\n",&N[z],&w[z]);
	}
	fclose(doping);
	for(z=0;z<NumLayers;z++){
	    N[z]=N[z]*1e6;
		w[z]=w[z]*1e-6;
	}
	G=new double[NumLayers];
	Gw=new double[NumLayers];
	Nwe=new double[NumLayers];
	Nw2e=new double[NumLayers];
	for(i=0;i<NumLayers;i++){
		G[i]=q*N[i]/die;
		G[0]=fabs(G[0]);		
		Gw[i]=G[i]*w[i];
		Nwe[i]=N[i]*w[i]/die;
		Nw2e[i]=N[i]*w[i]*w[i]/die;
	}
	
	Vtcheck=new double[NumLayers-1];
	Vtcheck[0]=0.5*(w[1]+Gw[1]/G[0])*Gw[1];
	for(i=1;i<NumLayers-1;i++){
		Vtcheck[i]=vtr(i);
	}
};
//profiler is the N_Layer electric field solver. Calculates field for voltagein as a
//N-I-P but is required as P-I-N so is then flipped and shifted so the first xposition is 0.
//PUBLIC
void device::profiler(double voltagein){
	double voltage;
	voltage=voltagein+Vbi;
	int i,numlayers;
	for(i=0;i<(NumLayers+1);i++){
		efield_x[i]=0;
		efield_e[i]=0;
	}
	double a,b,c,xlayer;
	if(voltage<Vtcheck[0] || NumLayers==2){
		//Two Layers
		numlayers=2;
		efield_x[0]=sqrt(2*voltage/(G[0]*(1+G[0]/G[1])));
		efield_x[1]=G[0]/G[1]*efield_x[0];
		efield_e[1]=G[0]*efield_x[0];
		efield_x[2]=efield_x[0]+efield_x[1];
		
	}
	else{
		//More than two layers
		for(i=1;i<(NumLayers-1);i++){
			if(voltage<Vtcheck[i] || NumLayers==i+2){
				numlayers=i+2;
				int k=i+2-1;//+2 converts to NumLayers depleted -1 to acount for change of base Matlab to C
				a=N[k]/die*(1-(N[k]/die)/(N[0]/die));
				b=2*N[k]/die*bmid(1, k-1);
				c=-(Wr(Nwe,1,k-1)*Wr(Nwe,1,k-1))/(N[0]/die)+(Wr(Nw2e,1,k-1))+cmid(k-1)-2*voltage/q;
				xlayer=(-b+sqrt(b*b-4*a*c))/(2*a);
				efield_x[0]=-(N[k]*xlayer/die+Wr(Nwe,1,k-1))/(N[0]/die);
				int j;
				for (j=0;j<k-1;j++){
					efield_e[j+1]=xlayer*G[k]+Wr(Gw,j+1,k-1);
				}
				efield_e[k]=xlayer*G[k];
				//xm=xlayer+Wr(w,1,k-1)+x[0];
				efield_x[2]=w[1];
				for (j=3;j<k+1;j++){
					efield_x[j]=Wr(w,1,j-1);
				}
				efield_x[k+1]=Wr(w,1,k-1)+xlayer;
				break;
			}
		}
	}
	efield_x[0]=-efield_x[0];
	double temp_efield_x[numlayers+1]={0};
	double temp_efield_e[numlayers+1]={0};
	i_max=numlayers;
	for(i=0;i<(numlayers+1);i++){
		temp_efield_x[i]=efield_x[numlayers-i];
		temp_efield_e[i]=efield_e[numlayers-i];
	}

	efield_x[0]=0;
	double temp;
	for(i=1;i<(numlayers+1);i++){
		temp=fabs(temp_efield_x[i-1]-temp_efield_x[i]);
		efield_x[i]=temp+efield_x[i-1];
	}
	
	for(i=0;i<(numlayers+1);i++){
		efield_e[i]=temp_efield_e[i];
	}

	width=efield_x[numlayers];
	//printf("%g\n",width);
};

