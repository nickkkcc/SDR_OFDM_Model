#include "DWT.h"
#include <C:\Users\Rise\source\repos\Model\wavelib-master\header\wavelib.h>
#include <math.h>


void DWT(Eigen::MatrixXd& Signal, const char* wname, Eigen::MatrixXd& cA, Eigen::MatrixXd& cD)
{
	int J = 1, N = Signal.size();
	wave_object obj = wave_init(wname);
	int cA_cD_length = floor((N + obj->hpd_len - 1) / 2);
	wt_object wt = wt_init(obj, "dwt", N, J);
	setDWTExtension(wt, "sym");
	setWTConv(wt, "direct");
	dwt(wt, &Signal(0));
	cA.resize(cA_cD_length, 1);
	cD.resize(cA_cD_length, 1);
	for (int i = 0,j = i + cA_cD_length; i < cA_cD_length; i++,j++)
	{
		cA(i, 0) = wt->output[i];
		cD(i, 0) = wt->output[j];
	}
	wave_free(obj);
	wt_free(wt);
};