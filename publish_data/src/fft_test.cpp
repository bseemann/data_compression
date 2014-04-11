#include <iostream>	// cout
#include <fstream>	// file operations
#include <cstdlib>	// malloc
#include <algorithm>	// sort
#include <cmath>
#include <complex>	// goes along with fftw
#include <fftw3.h>	

#define N 86400
#define Ncomplex 43201

using namespace std;

int main(){
	double *in_o, *in_p;
	fftw_complex *out;
	fftw_plan p;
	fftw_plan q;	
	
	// for calculating the threshold
	double *abs_fft;
	abs_fft = (double *)malloc(Ncomplex*sizeof(double));

	in_o = (double *)malloc(N*sizeof(double));
	in_p = (double *)malloc(N*sizeof(double));
	out = (fftw_complex *)fftw_malloc(Ncomplex*sizeof(fftw_complex));

	// construct plan
	p = fftw_plan_dft_r2c_1d(N, in_o, out,0);
	q = fftw_plan_dft_c2r_1d(N, out, in_p, 0);
	
	// read input data
	ifstream file ("office/processed/day1.txt",ifstream::in);
	if(file.is_open())
	{
		cout << "Reading file..." << endl;
		double x;
		for(int i=0;i<N;i++)
		{
			file >> x;
			(x>3500)? in_o[i]=1 : in_o[i]=0;
		//	cout << scientific << in[i] << endl;
		}
		file.close();
	}
	else cout << "Could not open file." << endl;

	// fft
	fftw_execute(p);

	ofstream ofs ("fft_transform.txt",ofstream::out | ofstream::trunc);
	ofs << "Output:" << endl;
	for(int i=0;i<Ncomplex;i++)
		//ofs << i << "\t" << creal(out[i]/N) << " " << cimag(out[i]/N) << endl;
		ofs << i << "\t" << out[i][0]/N << " " << out[i][1]/N << endl;
	ofs.close();

	ofs.open("abs_fft.txt",ofstream::out | ofstream::trunc);
	for(int i=0;i<Ncomplex;i++)
	{
		complex<double> x (out[i][0],out[i][1]);
		ofs << i << "\t" << abs(x) << endl;
		abs_fft[i] = abs(x);
	}
	ofs.close();

	ofs.open("abs_fft_sorted.txt",ofstream::out | ofstream::trunc);
	sort(abs_fft,abs_fft+(Ncomplex-1));
	for(int i=0;i<Ncomplex;i++)
		ofs << i << "\t" << abs_fft[i] << endl;
	ofs.close();
	
	ofs.open("fft_transform_thresholded.txt",ofstream::out | ofstream::trunc);
	double th = abs_fft[Ncomplex-16];
	for(int i=0;i<Ncomplex;i++)
	{
		complex<double> x (out[i][0],out[i][1]);
		if(abs(x) < th) out[i][0] = out[i][1] = 0;
		ofs << i << "\t" << out[i][0] << " " << out[i][1] << endl;
	}
	ofs.close();

	// ifft
	fftw_execute(q);

	ofs.open("reconstructed_signal.txt", ofstream::out | ofstream::trunc);
	for(int i=0;i<N;i++)
		ofs << i << "\t" << in_p[i]/N << endl;
	ofs.close();

	ofs.open("outlier.txt", ofstream::out | ofstream::trunc);
	for(int i=0;i<N;i++)
	{
		double x = (in_p[i]>0.5)?1:0;
		x = abs(x-in_o[i]);
		ofs << i << "\t" << x << endl;
	}
	ofs.close();


	fftw_destroy_plan(p); fftw_destroy_plan(q);
	fftw_free(out); 
	free(in_o); free(in_p); free(abs_fft);
	return 0;
}
