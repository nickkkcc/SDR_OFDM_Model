#include <iostream>
#include <Eigen/Dense>
#include <ctime>
#include <math.h>
#include <vector>
#include "ErrorCount.h"
#include "MessageClass.h"
#include "Randn.h"
#include <random>
#include <fftw3.h>
#include <complex>
#include <C:\Users\Rise\source\repos\Model\wavelib-master\header\wavelib.h>
#include "circ_shift.h"
#include "Hilbert.h"
#include <fstream>
#include "text_to_matrix.h"
#include "Net_Info.h"
#include <iomanip>



int main()
{
	Net_Info info;
	std::string p = "C:\\Users\\Rise\\Desktop\\НИР(VS)\\OFDMClassificationModel_Matlab\\";
	std::string PATHS_LW[3]{ p + "LW_layer_2.dat", p + "LW_layer_3.dat", p + "LW_layer_4.dat"};
	std::string PATH_IW = p + "IW.dat";
	std::string PATH_offset = p + "xoffset .dat";
	std::string PATH_gain = p + "gain .dat";
	std::string PATHS_b[4]{p + "bias_1.dat", p + "bias_2.dat" , p + "bias_3.dat", p + "bias_4.dat" };
	std::string PATH_ymin = p + "ymin .dat";
	std::string PATH_real = p + "realsignal.dat";
	std::string PATH_imag = p + "imagsignal.dat";
	std::string PATH_features = p + "features.dat";
	//getTextToMatrix(&PATH_ymin, info.ymin);
	//getTextToMatrix(PATH_gain, info.gain);
	//getTextToMatrix(PATH_IW, info.IW);
	//getTextToMatrix(PATH_offset, info.offset);

	//getTextToMatrix(PATH_features, info.features);
	//getTextToMatrix(PATHS_LW[0], info.LW2);
	//getTextToMatrix(PATHS_LW[1], info.LW3);
	//getTextToMatrix(PATHS_LW[2], info.LW4);
	//getTextToMatrix(PATHS_b[0], info.b1);
	//getTextToMatrix(PATHS_b[1], info.b2);
	//getTextToMatrix(PATHS_b[2], info.b3);
	//getTextToMatrix(PATHS_b[3], info.b4);
	

	



	srand(time(NULL));
	setlocale(LC_ALL, "RUS");
	double error = 0;
	int npack = 10, snr = 21, nsymb = 10, m = 1, nfft = 32;
	std::string duration = "ROWS";
	Eigen::MatrixXd features(22, 1);
	Eigen::MatrixXd gig(22, 1);
	for (size_t i = 0; i < snr + 1; i++)
	{
		for (size_t j = 1; j < npack + 1; j++)
		{
			
			Message A(nfft, nsymb, m);
            getTextToMatrix(PATH_real, PATH_imag, A.signal);

			
			A.message.resize(32, 1);
			for (int i = 0; i < A.signal.cols(); i++)
			{
				A.message = A.signal.col(i);

				A.getFeatures_14dwt3stat5cum(gig);
				features.col(i) = gig;
			}
			for (int i = 0; i < 10; i++)
			{
				std::cout << features(i,0) << " " << features(i,1) << " " << features(i,2) << " " << features(i,3) << std::endl;
			}
			A.Bi2Dec();
				//Модуляция:
			A.ModBPSK();
			
			//Показ модулированного сигнала

			A.Reshape(nfft, nsymb);

			/*A.ShowMessage();*/
			/*A.IFFTSHIFT(1);*/
			//Обратное преобразование Фурье
			A.DFT_inv("ROWS",nfft);
			A.getFeatures_14dwt3stat5cum(features);
			std::cout << features << std::endl;
			//A.ShowMessage();

			// АБГШ канал
			//A.AWGN(i + 10 * log10(2)); // Не использовать, не рабочий.

			//Прямое преобразование Фурье
			A.DFT_fwd("ROWS");
			/*A.ShowMessage();*/
			//Демодулирование
			/*A.DeModBPSK_EVK();*/
			A.DeModQPSK_EVK();
			//A.ShowMessage();
			//В двоичную форму
			
			A.Dec2Bi();
			A.ShowMessage();
		/*	error += ErrorCount(A, B);*/
		}
		std::cout << " Ошибка пакета при SNR = " << i << ": " << error / npack << std::endl;
	}

	system("pause");
};



