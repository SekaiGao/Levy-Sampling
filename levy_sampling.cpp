#pragma GCC optimize(3,"Ofast","inline")
#include<iostream>
#include<vector>
#include<ctime>
#include<cmath>
const long double pi=3.14159265358979323846;
using namespace std;
typedef long long ll;
typedef long double ld;

//erf inverse
inline ld erfInv(ld x)
{
   ld tt1, tt2, lnx, sgn;
   ld a=0.15449436008930206298828125;
   sgn=(x<0)?-1.:1.;
   
   x=(1-x)*(1+x);
   lnx=log(x);

   tt1=2./(pi*a)+0.5*lnx;
   tt2=lnx/a;

   return sgn*sqrtf(-tt1+sqrt(tt1*tt1-tt2));
}

//levy sampling
void levy_sampling(ld*u,int n,ld mu=0,ld c=1)
{
	int x0=time(NULL);
	ll a=16807,m=(1<<21);//LCG
	int i=0;
	while(i<n)
	{
		x0=(a*x0)%m;
		ld u0=ld(x0)/m;
		u[i++]=mu+c/(2.*pow(erfInv(u0),2));
	}
}


//caculate pdf
void pdf(ld*sample,ld*sta,int nsa,int nst,ld ls,ld rs)
{//sampling data, intervals, number of data, number of intervals, l/r bound
	ld dx=(rs-ls)/nst;
	for(int i=0;i<nsa;++i)
	{
		int index=(sample[i]-ls)/dx;
		if(index<nst)
		++sta[index];
	}
	ld dy=nsa*dx;
	for(int i=0;i<nst;++i)
	sta[i]/=dy;

}

int main()
{
	ld mu=0,c=1;
	int n=1000000;
	ld*levy_random_variable=new ld[n];
	//sampling
	levy_sampling(levy_random_variable,n,mu,c);


    int n_intervals = 1000;  
    ld *pdf_values = new ld[n_intervals](); 

    //pdf hist
    ld ls = 0; 
    ld rs = 100; 
    pdf(levy_random_variable, pdf_values, n, n_intervals, ls, rs);

	//write
	FILE*fp=fopen("levy.txt","w");
	for (int i = 0; i <n_intervals-1; ++i)
	fprintf(fp,"%Lf,",double(pdf_values[i]));
	fprintf(fp,"%Lf",double(pdf_values[n_intervals-1]));
	fclose(fp);


	delete[]levy_random_variable;
	delete[]pdf_values;

	return 0;
}