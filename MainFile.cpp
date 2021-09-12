//This Program will recognize the vowels using the provided test files or speech input

#include "stdafx.h"
#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#define PI 3.14159265

using namespace std;

long double ai[13];

long double ci[13];

long double R[13]={0};

long double F[11][6][13];

long double avgCi[6][13];

long double TokhuraDistance[6];

long double Ft[6][13];

long double Fr[6][13];

long double TokhuraWeights[13]={0, 1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

//finding DC shift and then correcting it
void DC_shiftCorrection(char *file)
{
	ifstream fin;
	ofstream fout;
	string line;
	
	char *file1="test.txt";
	
	long double amplitude,totalsum=0;
	long long sample_count=0;
	
	fin.open(file);
	
	while(!fin.eof())
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		totalsum+=amplitude;
		sample_count++;
	}	
	
	long double dc_shift;
	
	dc_shift=(double)totalsum/sample_count;
	
	fin.close();
	
	fin.open(file);
	fout.open(file1);
	
	while(!fin.eof())
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		amplitude=amplitude-dc_shift;
		fout<<amplitude<<endl;
	}
	fin.close();
	fout.close();
	
	fin.open(file1);
	fout.open(file);
	
	//copying the contents
	while(!fin.eof())
	{
		getline(fin,line);
		fout<<line<<endl;
	}
	
	fin.close();
	fout.close();
	
	if(remove(file1)!=0)
		perror("Error in removing file:");
}

//Normalization of the recordings upto max 10000 amplitude
void Normalization(char *file)
{
	ifstream fin;
	ofstream fout;
	string line;
	long double amplitude,max=INT_MIN;
	double norm_factor;

	fin.open(file);
	
	while(!fin.eof())
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		if(amplitude>max)
		max=amplitude;
	}
	fin.close();
	
	norm_factor=10000/max;
	
	fin.open(file);
	fout.open("Normalizedtest.txt");
	
	while(!fin.eof())
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		amplitude=amplitude*norm_factor;
		fout<<amplitude<<endl;
	}
	fin.close();
	fout.close();
	
	fin.open("Normalizedtest.txt");
	fout.open(file);
	
	//copying the contents of Normalizedtest file to test file
	while(!fin.eof())
	{
		getline(fin,line);
		fout<<line<<endl;
	}
	
	fin.close();
	fout.close();
	
	//remove the Normalizedtest file
	if(remove("Normalizedtest.txt")!=0)
		perror("Error removing in file:");
}

//Calculation of Ri values
void RiCalculation(char *file)
{
	ifstream fin;
	string line;
	long double buffer[400],amplitude;
	long double totalsum=0;

	fin.open(file);
	
	int sample_count=350;
		
	int i=0;	
	while(!fin.eof()&&sample_count--)
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		buffer[i]=amplitude;
		i++;
	}
		
	for(int i=0;i<=12;i++)
	{
		totalsum=0;
		for(int j=0;j<=320-1-i;j++)
		{
			totalsum+=buffer[j]*buffer[j+i];
		}
		R[i]=totalsum/320;
	}	
}

//Calculation of Ai values
void AiCalculation()
{
	int p=12;
	
	long double k[13],E[13];
	long double a[13][13];
	
	E[0]=R[0];
	k[1]=R[1]/E[0];
	a[1][1]=k[1];
	E[1]=(1-k[1]*k[1])*E[0];
	
	for(int i=2;i<=p;i++)
	{	
		long double sum=0;
		
	    for(int j=1;j<=i-1;j++)
	    {
	    	sum+=a[i-1][j]*R[i-j];
		}
		
		k[i]=(R[i]-sum)/E[i-1];
		
		a[i][i]=k[i];
		
		for(int j=1;j<=i-1;j++)
		{
			a[i][j]=a[i-1][j]-k[i]*a[i-1][i-j];
		}
		
		E[i]=(1-k[i]*k[i])*E[i-1];
	}
	
	for(int j=1;j<=p;j++)
	{	
		ai[j]=a[p][j];
	}
	
}

//Calculation of Ai coefficients
void CiCalculation()
{
	long double sum=0;
	
	ci[1]=ai[1];
	
	for(int m=2;m<=12;m++)
	{
		sum+=ai[m];
		
		for(int k=1;k<=m-1;k++)
		{	
			double x=(double)k/m;
			
			sum+=x*ci[k]*ai[m-k];		
		}	
		ci[m]=sum;
		sum=0;
	}
}

void ApplyRaisedSineWindow()
{
	for(int i=1;i<=12;i++)
	{
		long double w=(1+6*sin(PI*i/12));
		ci[i]=ci[i]*w;
	}
}

//Framing of training files this function will find 5 steady frames and will calculate Ci values
//for those frames
void Framing(char *file,long long fileno)
{
	ifstream fin;
	fin.open(file);	
	
	long double amplitude,total_sqr_amplitude=0,max=INT_MIN;
	long max_index;
	
	vector<long double> ste;
	int sample_count=0;
	string line;
	
	while(fin)
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		total_sqr_amplitude+=amplitude*amplitude;
		sample_count++;
		
		if(sample_count==320)
		{
			ste.push_back(total_sqr_amplitude/320);
			total_sqr_amplitude=0;
			sample_count=0;
		}
	}
	fin.close();
	
	for(unsigned int i=0;i<ste.size();i++)
	{
		if(max<ste[i])
		{
			max=ste[i];
			max_index=i;
		}
	}
	fin.open(file);
	
	long long skipline=320*(max_index-2);
	
	while(skipline--)
	{
		getline(fin,line);
	}
	long long frame_count=1;
	sample_count=0;
	ofstream fout;
	
	string str=to_string(frame_count)+".txt";
	const char *filename=&str[0];
	fout.open(filename);
	
	while(frame_count<=5)
	{
		if(sample_count==320)
		{
			frame_count++;
			fout.close();
			if(frame_count<=5)
			{
				string str=to_string(frame_count)+".txt";
				const char *filename=&str[0];
				
				fout.open(filename);
			}
			sample_count=0;
		}
		getline(fin,line);
		fout<<line<<endl;
		sample_count++;
	}
	
	for(long long i=1;i<=5;i++)
	{
		string str=to_string(i)+".txt";
		char *filename=&str[0];
		RiCalculation(filename);
		AiCalculation();
		CiCalculation();
		ApplyRaisedSineWindow();
		for(int j=1;j<=12;j++)
		{
			F[fileno][i][j]=ci[j];
		}
	}
	
	for(long long i=1;i<=5;i++)
	{
		string str=to_string(i)+".txt";
		const char *filename=&str[0];
		
		if(remove(filename)!=0)
		{
			perror("Error in removing files");
		}
	}
}

void AverageCi()
{
	for(int i=1;i<=5;i++)
	{
		for(int j=1;j<=12;j++)
		{	
		 long double sum=0;
		 
			for(int k=1;k<=10;k++)
			{
				sum+=F[k][i][j];
			}
			avgCi[i][j]=(long double)sum/10;
		}
	}
}

void creatingReferenceFile(char vowel)
{
	ofstream fout;
	string str=string("Reference file_")+vowel+".txt";	
	char *filename=&str[0];
	
	fout.open(filename);
	
	AverageCi();
	
	for(int i=1;i<=5;i++)
	{
		for(int j=1;j<=12;j++)
		{
			fout<<avgCi[i][j]<<" ";
		}
		fout<<endl;
	}
	fout.close();
}

void ReferenceFile_main()
{	
	cout<<"Creating Reference Files Please wait for 2 minutes...\n";

	char *vowels="aeiou";
	for(int k=0;k<5;k++)
	{	
		cout<<"\nCreating Reference File for vowel: "<<vowels[k]<<".......\n\n";	

		cout<<"Processing files: ";
		for(long long i=1;i<=10;i++)
		{	
			string str=string("Recordings//204101012_")+vowels[k]+"_"+to_string(i)+".txt";
			
			char *filename=&str[0];
		
		//	DC_shiftCorrection(filename);
	
		//	Normalization(filename);
	
			Framing(filename,i);

			cout<<i<<" ";
		}
	
		creatingReferenceFile(vowels[k]);

		cout<<"\n\nReference file created"<<endl;
	}
}

//Framing of test files this function will find 5 steady frames and will calculate Ci values for those
void FramingTestFiles(char *file)
{
	ifstream fin;
	fin.open(file);	
	
	long double amplitude,total_sqr_amplitude=0,max=INT_MIN;
	long max_index;
	
	vector<long double> ste;
	int sample_count=0;
	string line;
	
	while(fin)
	{
		getline(fin,line);
		stringstream s(line);
		s>>amplitude;
		total_sqr_amplitude+=amplitude*amplitude;
		sample_count++;
		
		if(sample_count==320)
		{
			ste.push_back(total_sqr_amplitude/320);
			total_sqr_amplitude=0;
			sample_count=0;
		}
	}
	fin.close();
	
	for(unsigned int i=0;i<ste.size();i++)
	{
		if(max<ste[i])
		{
			max=ste[i];
			max_index=i;
		}
	}
	fin.open(file);
	
	long long skipline=320*(max_index-2);
	
	while(skipline--)
	{
		getline(fin,line);
	}
	long long frame_count=1;
	sample_count=0;
	ofstream fout;
	
	string str=to_string(frame_count)+".txt";
	const char *filename=&str[0];
	fout.open(filename);
	
	while(frame_count<=5)
	{
		if(sample_count==320)
		{
			frame_count++;
			fout.close();
			if(frame_count<=5)
			{
				string str=to_string(frame_count)+".txt";
				const char *filename=&str[0];
				
				fout.open(filename);
			}
			sample_count=0;
		}
		getline(fin,line);
		fout<<line<<endl;
		sample_count++;
	}
	
	for(long long i=1;i<=5;i++)
	{
		string str=to_string(i)+".txt";
		char *filename=&str[0];
		RiCalculation(filename);
		AiCalculation();
		CiCalculation();
		ApplyRaisedSineWindow();
		for(int j=1;j<=12;j++)
		{
			Ft[i][j]=ci[j];
		}
	}
	
	for(long long i=1;i<=5;i++)
	{
		string str=to_string(i)+".txt";
		const char *filename=&str[0];
		
		if(remove(filename)!=0)
		{
			perror("Error in removing files");
		}
	}
}

void CalculateTokhuraDistance()
{
	ifstream fin;

	string word;
	long double CiValues;
	char *vowels="aeiou";
	
	for(int k=0;k<5;k++)
	{
	string str=string("Reference file_")+vowels[k]+".txt";	
	char *filename=&str[0];
	fin.open(filename);
	
	for(int i=1;i<=5;i++)
	{
		for(int j=1;j<=12;j++)
		{
			fin>>word;
			stringstream s(word);
			s>>CiValues;
			
			Fr[i][j]=CiValues;
		}
	}
	
	fin.close();
	
	long double avgTokhura=0;
	for(int i=1;i<=5;i++)
	{
		long double sum=0;
		for(int j=1;j<=12;j++)
		{
			long double diff=Ft[i][j]-Fr[i][j];
			sum+=TokhuraWeights[j]*diff*diff;
		}
		avgTokhura+=sum;
	}
	avgTokhura=(double)avgTokhura/5;
	TokhuraDistance[k+1]=avgTokhura;
	}
	
}

int main()
{	
	ReferenceFile_main();

	int count=0,choice;
	char *vowels="aeiou";
	cout<<"\nRecognizing vowels in the test files: \n\n";

	cout<<"Choose from the following options:";
	cout<<"\n\n1.Take input from the test files already recorded\n";
	cout<<"\n2.Take input from the microphone of device(Recording Module Implementation)\n";
	cout<<"\nEnter your choice(1/2): ";
	cin>>choice;
	cout<<endl;

	if(choice==1)
	{
	for(int k=0;k<5;k++)
	{
		for(long long x=11;x<=20;x++)
		{
			string str=string("Recordings//204101012_")+vowels[k]+"_"+to_string(x)+".txt";
			cout<<"Real Vowel is "<<vowels[k]<<"; ";
		long double min=INT_MAX;
		int min_index;
		
		char *filename=&str[0];
		
		DC_shiftCorrection(filename);
	
		Normalization(filename);
	
		FramingTestFiles(filename);
		
		CalculateTokhuraDistance();
		
		for(int i=1;i<=5;i++)
		{
			if(TokhuraDistance[i]<min)
			{
				min=TokhuraDistance[i];
				min_index=i;
			}
			
		}
		
		switch(min_index)
		{
			case 1: cout<<"Vowel Detected is a"<<endl;
			break;
			case 2: cout<<"Vowel Detected is e"<<endl;
			break;
			case 3: cout<<"Vowel Detected is i"<<endl;
			break;
			case 4: cout<<"Vowel Detected is o"<<endl;
			break;
			case 5: cout<<"Vowel Detected is u"<<endl;
		}
		
		if(k==min_index-1) count++;
	}
}
	cout<<"Accuracy: "<<count<<"/100\n";
	}
	else
	{
		system("Recording_Module.exe 3 MySample.wav MySample.txt");

		string str="MySample.txt";

		long double min=INT_MAX;
		int min_index;
		
		char *filename=&str[0];
		
		DC_shiftCorrection(filename);
	
		Normalization(filename);
	
		FramingTestFiles(filename);
		
		CalculateTokhuraDistance();
		
		for(int i=1;i<=5;i++)
		{
			if(TokhuraDistance[i]<min)
			{
				min=TokhuraDistance[i];
				min_index=i;
			}
			
		}
		
		switch(min_index)
		{
			case 1: cout<<"Vowel Detected is a"<<endl;
			break;
			case 2: cout<<"Vowel Detected is e"<<endl;
			break;
			case 3: cout<<"Vowel Detected is i"<<endl;
			break;
			case 4: cout<<"Vowel Detected is o"<<endl;
			break;
			case 5: cout<<"Vowel Detected is u"<<endl;
		}

	}
	return 0;
}