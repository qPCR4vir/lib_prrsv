//=============================================================================
// Module:        talign_test.cpp
// Project:       Diploma Thesis - Probe Selection for DNA Microarrays
// Type:          implementation - Thermodynamic Alignment.
// Language:      c++
// Compiler:      microsoft visual c++ 6.0, unix/linux gcc
// System/OS:     Windows 32, Sun solaris, Linux, other unix systems (untested)
// Database:      none
// Description:   ThermAlign - Thermodynamic Alignment Test Program
// Author:        kaderali
// Date:          9/2000
// Copyright:     (c) L. Kaderali, 9/2000
//
// Revision History
// $              00sep04 : created LK
// #$
//=============================================================================

//#pragma pack(1)

using namespace std;

#include "stdafx.h"
#include <string>
#include <stdlib.h>
#include <iostream>
#include "nnparams.h"
#include "thermalign.h"
#include "dinkelbach.h"
#include <stdio.h>

using namespace std;

#define cfilter 7.6

/*
int main(int argc, char* argv[])
{
	printf("TermAlign - Alignment of two Sequences using nearest neighbor"
 	         "thermodynamics\n            as weight function.\n");
	printf("=============================================================="
		   "=================\n"); 	
	char seq1[20000];
	char seq2[20000];

        int eins = 0;
        int zwei = 0;
        int drei = 0;
        int vier = 0;
        int fuenf = 0;
        int sechs = 0;
        int sieben = 0;
        int acht = 0;
        int neun = 0;
        int zehn = 0;
        int gzehn = 0;

	int seq_num, seq_length;
	cout << "How many sequences? ";
	cin >> seq_num;
	cout << "Determine sequence length: ";
	cin >> seq_length;

        srand(static_cast<unsigned>(time(0)));
        int sequence1;

 	for(int e = 0; e < seq_num; e++)
	{
                //generation of 2 random sequences with a length of 15 nucleotides
                seq1[0] = '$';
                seq2[0] = '$';

                for(int i=1; i < seq_length+1; i++)
                {
                  sequence1 = rand() % 4;
                  if(sequence1 == 0) seq1[i]='A';
                  if(sequence1 == 1) seq1[i]='T';
                  if(sequence1 == 2) seq1[i]='G';
                  if(sequence1 == 3) seq1[i]='C';
                }

                for(int h=1; h < seq_length+1; h++)
                {
                  sequence1 = rand() % 4;
                  if(sequence1 == 0) seq2[h]='A';
                  if(sequence1 == 1) seq2[h]='T';
                  if(sequence1 == 2) seq2[h]='G';
                  if(sequence1 == 3) seq2[h]='C';
                }

                seq1[seq_length+1] = '$';
                seq2[seq_length+1] = '$';
                seq1[seq_length+2] = '\0';
                seq2[seq_length+2] = '\0';

		float conc = 0.000001;

		PNNParams myParams = new CNNParams();
	        myParams->InitParams(conc,1);

                GGAlign mydGAlign = new GAlign(strlen(seq1),strlen(seq2),myParams); //verändert von ML!
                mydGAlign->InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));

                dinkelbach myDinkel(myParams, mydGAlign); //verändert von ML!

		CThermAlign myAlignment(strlen(seq1),strlen(seq2),myParams);
                myAlignment.InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));
                myAlignment.InitBorder(); 	
                myAlignment.CalculateTable();
                //std::cout<<endl<<"Sequence 1: "<<seq1<<endl<<"Sequence 2: "<<seq2 <<endl<<endl<<flush;

                myAlignment.OutputLocalAlignment(std::cout);

                float TempK = myAlignment.GetMeltingTempK(myAlignment.maxloci,myAlignment.maxlocj); //von ML
                int num_iteration = myDinkel.iteration(TempK);    //von ML
                //cout << "Anzahl der nötigen Iterationen: " <<num_iteration<<endl;
                if(num_iteration == 1) eins++;
                if(num_iteration == 2) zwei++;
                if(num_iteration == 3) drei++;
                if(num_iteration == 4) vier++;
                if(num_iteration == 5) fuenf++;
                if(num_iteration == 6) sechs++;
                if(num_iteration == 7) sieben++;
                if(num_iteration == 8) acht++;
                if(num_iteration == 9) neun++;
                if(num_iteration == 10) zehn++;
                if(num_iteration > 10) gzehn++;

		//std::cout<<endl<<flush;
	
	}
        cout << "1  = " << eins << endl;
        cout << "2  = " << zwei << endl;
        cout << "3  = " << drei << endl;
        cout << "4  = " << vier << endl;
        cout << "5  = " << fuenf << endl;
        cout << "6  = " << sechs << endl;
        cout << "7  = " << sieben << endl;
        cout << "8  = " << acht << endl;
        cout << "9  = " << neun << endl;
        cout << "10 = " << zehn << endl;
        cout << "g10= " << gzehn << endl;

	return 0;
}
*/

int main(int argc, char* argv[])
{
	printf("TermAlign - Alignment of two Sequences using nearest neighbor"
 	         "thermodynamics\n            as weight function.\n");
	printf("=============================================================="
		   "=================\n");
	char seq1[3000];
	char seq2[3000];

	bool end = false;

	std::cout<<"q to abort!"<<endl;

	while (!end)
	{
		std::cout<<endl<<"----"<<endl<<"Sequence 1: ";
		cin>>seq1;
		if (seq1[0]=='q')
			break;
		cout<<"Sequence 2: ";
		std::cin>>seq2;
		float conc1,conc2;
		float k;
		std::cout<<"Concentration 1 (e.g. 0.000001): ";
		std::cin>>conc1;
		std::cout<<"Concentration 2 (e.g. 0.000001): ";
		std::cin>>conc2;
		//conc2=0;
		
		std::cout<<"Total salt concentration  (e.g. 1): ";
		std::cin>>k;

		PNNParams myParams = new CNNParams();
	        myParams->InitParams(conc1,conc2,k,SALT_METHOD_SANTALUCIA);

                GGAlign mydGAlign = new GAlign(strlen(seq1),strlen(seq2),myParams);
                mydGAlign->InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));

                dinkelbach myDinkel(myParams, mydGAlign);

		CThermAlign myAlignment(strlen(seq1),strlen(seq2),myParams);
                myAlignment.InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));
                myAlignment.InitBorder();
                myAlignment.CalculateTable();
                //std::cout<<endl<<"Sequence 1: "<<seq1<<endl<<"Sequence 2: "<<seq2 <<endl<<endl<<flush;

                //myAlignment.OutputLocalAlignment(std::cout);

                double TempK = myAlignment.GetMeltingTempK(myAlignment.maxloci,myAlignment.maxlocj);
                //cout << "Temp before: " <<TempK<<endl;
		
		int num_iteration = myDinkel.iteration(TempK);    //von ML
                cout << "number of iterations: " <<num_iteration<<endl;
                cout << "Temperature: " <<mydGAlign->GetMeltingTempC(mydGAlign->maxloci,mydGAlign->maxlocj)<<endl;

		std::cout<<endl<<flush;

	}

	return 0;
}

/*
char getComp(char base) {
	if(base == 'A') return 'T';
	if(base == 'C') return 'G'; 
	if(base == 'G') return 'C'; 
	if(base == 'T') return 'A';
	return ' '; 
}
void getComplement(char * cadIni, char * cadOut)
{
	int len = strlen(cadIni);
	for(int i=0;i<len;i++) {
		cadOut[i] = getComp(cadIni[i]);
	}
	cadOut[len] = 0;
}
int countGCBases(char * seq ) {
	int len = strlen(seq);
	int count = 0;
	for(int i=0;i<len;i++) {
		if(seq[i] == 'C' || seq[i] == 'G') {
			count++;
		}
	}
	return count;
}

int countMismatches(char * seq1, char * seq2) {
	int answer = 0;
	for(int i=0;i<strlen(seq1);i++) {
		if(seq1[i] != getComp(seq2[i])) {
			answer++;
		}
	}
	return answer;
}

int calculateTemperature(char * sequen1, char * sequen2, float conc, float k,int method, PNNParams myParams, float * temp, float * gcContent)
{
	*gcContent = countGCBases(sequen1)+countGCBases(sequen2);
	float totalLength = strlen(sequen1)+strlen(sequen2);
	*gcContent /= totalLength;
	char seq1[255];
	char seq2[255];
	strcpy(seq1,"$");
	strcat(seq1,sequen1);
	strcat(seq1,"$");
	strcpy(seq2,"$");
	strcat(seq2,sequen2);
	strcat(seq2,"$");
	float singleConc = conc/2;
	myParams->InitParams(singleConc,singleConc,k,method);
	
	CThermAlign myAlignment(strlen(seq1),strlen(seq2),myParams);
	myAlignment.InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));
	myAlignment.InitBorder();
	myAlignment.CalculateTable();
 

	double TempK = myAlignment.GetMeltingTempK(myAlignment.maxloci,myAlignment.maxlocj);
	GAlign mydGAlign (strlen(seq1),strlen(seq2),myParams);
	mydGAlign.InitStrings(seq1,seq2,strlen(seq1),strlen(seq2));

	dinkelbach myDinkel(myParams, &mydGAlign);

        int num_iteration = myDinkel.iteration(TempK);	
	*temp = mydGAlign.GetMeltingTempC(mydGAlign.maxloci,mydGAlign.maxlocj);
	return num_iteration;
}
int main(int argc, char* argv[])
{
	char seq1[255];
	char seq2[255];
	char line[1000],*tok;
	int id;
	float gcContent;
	float conc;
	float k;
	float tmExp;
	float temp;
	float temp2;
	FILE * fin, *fout;
	fin = fopen("input.csv","rt");
	fout = fopen("output.csv","wt");
	PNNParams myParams = new CNNParams ();
		fgets(line,1000,fin);
	while (!feof(fin))
	{
		tok = strtok(line,",");
		id = atoi(tok);
		tok = strtok(NULL,",");
		strcpy(seq1,tok);
		tok = strtok(NULL,",");
		strcpy(seq2,tok);
		if(strcmp(seq2,"COMPLEMENT")==0) {
			getComplement(seq1,seq2);
		}
		
		tok = strtok(NULL,",");
		gcContent = atof(tok);
		tok = strtok(NULL,",");
		conc = atof(tok);
		tok = strtok(NULL,",");
		k = atof(tok);
		tok = strtok(NULL,",");
		tmExp = atof(tok);
		int num_iteration = calculateTemperature(seq1, seq2, conc, k, SALT_METHOD_SANTALUCIA,  myParams, &temp,&gcContent);
		calculateTemperature(seq1, seq2, conc, k, SALT_METHOD_OWCZARZY, myParams, &temp2,&gcContent);
		int mismatches = countMismatches(seq1,seq2);
		fprintf(fout,"%d,%s,%s,%d,%.2f,%d,%.3f,%.3f,%d\n",id,seq1,seq2,(strlen(seq1)+strlen(seq2))/2,100*gcContent,mismatches,temp,temp2,num_iteration);
		fgets(line,1000,fin);
	}
	fclose(fin);
	fclose(fout);
	return 0;
}
*/
