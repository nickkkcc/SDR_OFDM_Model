#pragma once
#include<iostream>
#include<Eigen/Dense>
#include<vector>
#include <unsupported/Eigen/FFT>
#include <fftw3.h>

class Message
{
public:
	int nfft; // Окно преобразования Фурье.
	int nsymb; // Кол-во символов
	int m;  // Кол-во бит на символ
	int rows; //Кол-во строк 
	int cols;
	Eigen::MatrixXcd message; //Экземпляр матрицы Eigen.
	Eigen::MatrixXcd signal;
	fftw_plan plan_fwd;
	fftw_plan plan_inv;


	Message();
	Message(int, int, int);
	void ShowMessage();
	Eigen::MatrixXcd& Bi2Dec();
	Eigen::MatrixXcd& Dec2Bi();
	Eigen::MatrixXcd& ModBPSK();
	Eigen::MatrixXcd& DeModBPSK();
	Eigen::MatrixXcd& ModQPSK();
	Eigen::MatrixXcd& DeModQPSK();
	Eigen::MatrixXcd& Mod16QAM();
	Eigen::MatrixXcd& DeMod16QAM();
	Eigen::MatrixXcd& Reshape(int, int);
	Eigen::MatrixXcd& DFT_inv(std::string, int = 1);
	Eigen::MatrixXcd& DFT_fwd(std::string, int = 1);
	Eigen::MatrixXcd& DeModBPSK_EVK();
	Eigen::MatrixXcd& DeModQPSK_EVK();
	Eigen::MatrixXcd& DeMod16QAM_EVK();
	//Eigen::MatrixXcd& AWGN(int);
	int Size();
	int Rows();
	int Cols();
	double GetReal(int i, int j);
	double GetImag(int i, int j);
	Eigen::MatrixXcd& IFFTSHIFT(int dir);
	Eigen::MatrixXcd& IFFTSHIFT();
	Eigen::MatrixXcd& FFTSHIFT(int dir);
	Eigen::MatrixXcd& FFTSHIFT();
	void getFeatures_14dwt3stat5cum(Eigen::MatrixXd& features);

};
double ErrorCount(Message&, Message&);