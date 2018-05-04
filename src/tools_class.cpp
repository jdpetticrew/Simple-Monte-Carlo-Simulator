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
   tools_class.cpp contains the class implimentation of the tools class for the SMC
   The tools class calculates the scattering rates and probabilities for the SMC

   tools.h contains the class definition

   Jonathan Petticrew, University of Sheffield, 2017.
 */

#include "tools.h"
#include <conio.h>
#include <math.h>
#include <stdio.h>
//Constructs and zeros the probability arrays
tools::tools(SMC *input) : constants(input){
	int i,j;
	for (i=0; i<PB_DIM1_SZ; i++)
		for(j=0; j<PB_DIM2_SZ; j++)
		{
			pb[i][j]=0;
			pb2[i][j]=0;
		}
};
//Calculates the Scattering Probabilities
int tools::scattering_probability(){ //calculates the scattering probabilities contained in pb[][] and pb2[][]
	FILE *fp_rate, *fp_pb;
	bool ratefail=false;
	bool pbfail=false;
	unsigned int j,inputkey=0,GoAhead=1;
	e_rate();
	h_rate();
	int x = constants->Get_NUMPOINTS();
	rtotal=pb[0][x]+pb[1][x]+pb[2][x]; //rtotal is the total interaction rate for electrons at the Max energy declaired in class SMC
	rtotal2=pb2[0][x]+pb2[1][x]+pb2[2][x]; //rtotal 2 is similar for holes

	if ((fp_rate=fopen("scattering_rates.txt","w"))==NULL)
	{    printf("Cannot oen file \"scattering_rates.txt\"\n");
	     ratefail=true;}
	if ((fp_pb=fopen("scattering_pb.txt","w"))==NULL)
	{    printf("Cannot oen file \"scattering_pb.txt\"\n");
	     pbfail=true;}
	if(ratefail||pbfail)
	{    printf("Output files could not be opened. Would you like to continue (y/n)\n");
	     do
	     {inputkey=_getch();}
	     while(inputkey!='y' && inputkey!='Y' && inputkey!='n' && inputkey!='N');
	     if (inputkey=='y' || inputkey=='Y')
		     GoAhead=2;
	     else GoAhead=0; }
	/****CHANGES THE RATES INTO PROBABILITIES****/
	if (GoAhead)
	{    for(j=0; j<=x; j++)
	     {
		     if(!ratefail)
			     fprintf(fp_rate,"%f, %e, %e, %e, %e, %e, %e\n",j*0.001,pb[0][j],pb[1][j],pb[2][j],pb2[0][j],pb2[1][j],pb2[2][j]);
		     pb[0][j]= pb[0][j]/rtotal;
		     pb[1][j]=pb[0][j]+pb[1][j]/rtotal;
		     pb[2][j]=pb[1][j]+pb[2][j]/rtotal;
		     pb2[0][j]= pb2[0][j]/rtotal2;
		     pb2[1][j]=pb2[0][j]+pb2[1][j]/rtotal2;
		     pb2[2][j]=pb2[1][j]+pb2[2][j]/rtotal2;
		     if(!pbfail)
			     fprintf(fp_pb,"j=%d, %e, %e, %e, %e, %e, %e\n",j,pb[0][j],pb[1][j],pb[2][j],pb2[0][j],pb2[1][j],pb2[2][j]);
	     }}
	if(!ratefail) fclose(fp_rate);
	if(!pbfail) !fclose(fp_pb);
	return(GoAhead);
	// returns 0 for no output and user wants to quit
	//returns 1 if everything ok
	//returns 2 if no output but user wants to continue
};

//calculates the electron scattering rates
void tools::e_rate(){
	unsigned int i;
	double n;
	double Egap;
	double x,x1;
	double e_para, e_para2;
	e_para=constants->Get_N()/(constants->Get_e_meanpath()*(2*constants->Get_N()+1));
	e_para2=(constants->Get_N()+1)/(constants->Get_e_meanpath()*(2*constants->Get_N()+1));
	int j = constants->Get_NUMPOINTS();

	for (i=0; i<=j; i++)
	{    n=i*constants->Get_q()*0.001;
	     Egap=(n-constants->Get_e_Eth())/(constants->Get_e_Eth());
	     pb[0][i]=e_para*sqrt((2*(n+constants->Get_hw()))/constants->Get_e_mass());//rate of phonon absorption

	     if (n>constants->Get_hw())
		     pb[1][i]=e_para2*sqrt((2*(n-constants->Get_hw()))/constants->Get_e_mass()); //rate of phonon emission
	     else pb[1][i]=0;

	     if (Egap>0)
	     { x=my_pow(Egap,constants->Get_e_gamma());
	       x1=constants->Get_e_Cii();
	       pb[2][i]=x*constants->Get_e_Cii();         //rate of impact ionization
	     }
	     else pb[2][i]=0; }
};
//calculates the hole scattering rates
void tools::h_rate(){
	unsigned int i;
	double n;
	double Egap;
	double h_para, h_para2;
	double x2,x1;
	h_para=constants->Get_N()/(constants->Get_h_meanpath()*(2*constants->Get_N()+1));
	h_para2=(constants->Get_N()+1)/(constants->Get_h_meanpath()*(2*constants->Get_N()+1));
	int x=constants->Get_NUMPOINTS();
	for (i=0; i<=x; i++)
	{    n=i*constants->Get_q()*0.001;
	     Egap=(n-constants->Get_h_Eth())/(constants->Get_h_Eth());
	     pb2[0][i]=h_para*sqrt((2*(n+constants->Get_hw()))/constants->Get_h_mass());//rate of phonon absorption

	     if (n>constants->Get_hw())
		     pb2[1][i]=h_para2*sqrt((2*(n-constants->Get_hw()))/constants->Get_h_mass()); //rate of phonon emission
	     else pb2[1][i]=0;

	     if (Egap>0)
	     { x2=my_pow(Egap,constants->Get_h_gamma());
	       x1=constants->Get_h_Cii();
	       pb2[2][i]=x1*x2;         //rate of impact ionization
	     }
	     else pb2[2][i]=0; }
};

//Get scattering rate of electrons total
double tools::Get_rtotal(){
	return rtotal;
};
//Pb get function electrons
double tools::Get_pb(int i, int j){
	return pb[i][j];
};
//Get scattering rate of holes total
double tools::Get_rtotal2(){
	return rtotal2;
};
//Pb get function holes
double tools::Get_pb2(int i, int j){
	return pb2[i][j];
};
/*my_pow is used to fix a bug with the pow function in my compiler. My compiler TDM-GCC 4.9.2 has an over
   accuracy problem with the pow function, if the value is not stored straight away it might return a
   different value.*/
double tools::my_pow(double base, double exponent)
{
	double output=pow(base,exponent);
	return output;
};
