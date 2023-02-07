#include "Hilbert.h"
#include <math.h>
#include <fftw3.h>

void Hilbert(Eigen::MatrixXcd& Input, Eigen::MatrixXcd& Output)
{
	Eigen::MatrixXd Temp(Input.rows(), Input.cols());
	Temp = Input.real();
	int N = Input.size();
	int fp_start = 1, fp_end = (floor(N / 2) + (N % 2)) - 1;
	int fn_start = ceil(N / 2) + abs(N % 2 - 1), fn_end = N - 1;
	fftw_plan plan_fft = fftw_plan_dft_1d(N, (fftw_complex*)& Input(0), (fftw_complex*)& Output(0), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan_fft);
	fftw_destroy_plan(plan_fft);
	Eigen::MatrixXcd temp(Input.rows(),Input.cols());
	for (int i = 0; i < temp.size(); i++)
	{
		temp(i) = { 0,0 };
	}
	for (int i = fp_start; i <= fp_end; i++)
	{
		temp(i) = Output(i) + Output(i);

	}
	plan_fft = fftw_plan_dft_1d(N, (fftw_complex*)& temp(0), (fftw_complex*)& Output(0), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_fft);
	fftw_destroy_plan(plan_fft);
	for (int i = 0; i < N; i++)
	{
		Output(i) = { Temp(i), Output(i).imag() / N };

	}
};