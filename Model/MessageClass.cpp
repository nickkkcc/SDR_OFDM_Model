#include<iostream>
#include "MessageClass.h"
#include <math.h>
#include <vector>
#include "Randn.h"
#include "circ_shift.h"
#include <C:\Users\Rise\source\repos\Model\wavelib-master\header\wavelib.h>
#include "DWT.h"
#include "HIlbert.h"


Message::Message() //Конструктор по умолчанию.
{
	std::cout << "Default modulation is BPSK (nffr = 32, nsymb = 10, m = 2)" << std::endl;
	nsymb = 10;
	nfft = 32;
	m = 2;
	int value = nfft * nsymb;
	message = Eigen::MatrixXcd().real().Random(value, m);
	rows = message.rows();
	cols = message.cols();

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (message.real()(i, j) > 0)
			{
				message(i, j) = 1;
			}
			else
			{
				message(i, j) = 0;
			}

		}
	}


};
Message::Message(int nfft_in, int nsymb_in, int m_in) //Конструктор
{
	nfft = nfft_in;
	nsymb = nsymb_in;
	m = m_in;
	int value = nfft * nsymb;
	message = Eigen::MatrixXcd().real().Random(value, m);
	rows = message.rows();
	cols = message.cols();

	for (int i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < cols; j++)
		{
			if (message.real()(i, j) > 0)
			{
				message.real()(i, j) = 1;
			}
			else
			{
				message.real()(i, j) = 0;
			}

		}
	}



};
Eigen::MatrixXcd& Message::Bi2Dec() //Перево в десятичные числа, исходя из вида модуляции
{
	
	switch (m)
	{
	case 2:
	{
		for (int i = 0; i < rows; i++)
		{
			message.real()(i, 0) = message.real()(i, 1) +
				message.real()(i, 0) * 2;
		}
		message.conservativeResize(Eigen::NoChange_t::NoChange, 1);
		cols = message.cols();
		break;
	}
	case 4:
	{
		for (int i = 0; i < rows; i++)
		{
			message.real()(i, 0) = message.real()(i, 3) +
				message.real()(i, 2) * 2 +
				message.real()(i, 1) * 4 +
				message.real()(i, 0) * 8;
		}
		message.conservativeResize(Eigen::NoChange_t::NoChange, 1);
		cols = message.cols();
		break;
	}
	case 6:
	{
		for (int i = 0; i < rows; i++)
		{
			message.real()(i, 0) = message.real()(i, 5) +
				message.real()(i, 4) * 2 +
				message.real()(i, 3) * 4 +
				message.real()(i, 2) * 8 +
				message.real()(i, 1) * 16 +
				message.real()(i, 0) * 32;
		}
		message.conservativeResize(Eigen::NoChange_t::NoChange, 1);
		cols = message.cols();
		break;
	}
	case 8:
	{
		for (int i = 0; i < rows; i++)
		{
			message.real()(i, 0) = message.real()(i, 7) +
				message.real()(i, 6) * 2 +
				message.real()(i, 5) * 4 +
				message.real()(i, 4) * 8 +
				message.real()(i, 3) * 16 +
				message.real()(i, 2) * 32 +
				message.real()(i, 1) * 64 +
				message.real()(i, 0) * 128;
		}
		message.conservativeResize(Eigen::NoChange_t::NoChange, 1);
		cols = message.cols();
		break;
	}
	}
	return message;
}
Eigen::MatrixXcd& Message::Dec2Bi()//Перевод двоичное представление
{
	
	int value = 0;
	double temp = rows * cols;
	message = message.reshaped<double, double>(temp, 1);
	message.conservativeResize(Eigen::NoChange_t::NoChange, m);
	for (int i = 0; i < rows * cols; i++)
	{
		value = message.real()(i, 0);
		for (int j = m - 1; j > -1; j--)
		{
			message(i, j) = value % 2;
			value = (value - message.real()(i, j)) / 2;
		}
	};
	cols = m;
	return message;
}
Eigen::MatrixXcd& Message::ModBPSK() // Модуляция BPSK
{

	for (int i = 0; i < rows; i++)
	{
		if (message.real()(i, 0) == 1)
		{
			message.real()(i, 0) = -1;
		}
		if (message.real()(i, 0) == 0)
		{
			message.real()(i, 0) = 1;
		}
	}
	return message;
}
Eigen::MatrixXcd& Message::DeModBPSK() //Демодуляция BPSK
{
	int size = rows * cols;
	message = message.reshaped<double, double>(size, 1);
	cols = message.cols();
	rows = message.rows();
	for (int i = 0; i < size; i++)
	{
		if (message.real()(i, 0) >= 0)
		{
			message.real()(i, 0) = 0;
			message.imag()(i, 0) = 0;
		}
		else
		{
			message.real()(i, 0) = 1;
			message.imag()(i, 0) = 0;
		}
	}
	return message;
}
Eigen::MatrixXcd& Message::ModQPSK() //Модуляция QPSK
{
	
	double value = 0.70710548251124;
	for (int i = 0; i < rows; i++)
	{
		int temp = message.real()(i, 0);
		switch (temp)
		{
		case 0:
		{
			message.real()(i, 0) = value;
			message.imag()(i, 0) = value;
			break;
		}
		case 1:
		{
			message.real()(i, 0) = value;
			message.imag()(i, 0) = -1 * value;
			break;
		}
		case 2:
		{
			message.real()(i, 0) = -1 * value;
			message.imag()(i, 0) = value;
			break;
		}
		case 3:
		{
			message.real()(i, 0) = -1 * value;
			message.imag()(i, 0) = -1 * value;
			break;
		}
		}
	}
	return message;
}
Eigen::MatrixXcd& Message::DeModQPSK()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (message.real()(i, j) >= 0 && message.imag()(i, j) >= 0)
			{
				message.real()(i, j) = 0;
				message.imag()(i, j) = 0;
				continue;
			}
			if (message.real()(i, j) >= 0 && message.imag()(i, j) < 0)
			{
				message.real()(i, j) = 1;
				message.imag()(i, j) = 0;
				continue;
			}
			if (message.real()(i, j) < 0 && message.imag()(i, j) >= 0)
			{
				message.real()(i, j) = 2;
				message.imag()(i, j) = 0;
				continue;
			}
			if (message.real()(i, j) < 0 && message.imag()(i, j) < 0)
			{
				message.real()(i, j) = 3;
				message.imag()(i, j) = 0;
				continue;
			}
		}

	}
	return message;
}
Eigen::MatrixXcd& Message::Mod16QAM()
{
	message = message.reshaped<double, double>(rows * cols, 1);
	cols = message.cols();
	rows = message.rows();
	double temp = rows * cols;
	for (int i = 0; i < temp; i++)
	{
		//std::cout << "Iteration number: " << i << std::endl;
		int temp = message.real()(i, 0);
		switch (temp)
		{
		case 0:
		{
			message.real()(i, 0) = -3;
			message.imag()(i, 0) = -3;
			break;
		}
		case 1:
		{
			message.real()(i, 0) = -3;
			message.imag()(i, 0) = -1;
			break;
		}
		case 2:
		{
			message.real()(i, 0) = -3;
			message.imag()(i, 0) = 3;
			break;
		}
		case 3:
		{
			message.real()(i, 0) = -3;
			message.imag()(i, 0) = 1;
			break;
		}
		case 4:
		{
			message.real()(i, 0) = -1;
			message.imag()(i, 0) = -3;
			break;
		}
		case 5:
		{
			message.real()(i, 0) = -1;
			message.imag()(i, 0) = -1;
			break;
		}
		case 6:
		{
			message.real()(i, 0) = -1;
			message.imag()(i, 0) = 3;
			break;
		}
		case 7:
		{
			message.real()(i, 0) = -1;
			message.imag()(i, 0) = 1;
			break;
		}
		case 8:
		{
			message.real()(i, 0) = 3;
			message.imag()(i, 0) = -3;
			break;
		}
		case 9:
		{
			message.real()(i, 0) = 3;
			message.imag()(i, 0) = -1;
			break;
		}
		case 10:
		{
			message.real()(i, 0) = 3;
			message.imag()(i, 0) = 3;
			break;
		}
		case 11:
		{
			message.real()(i, 0) = 3;
			message.imag()(i, 0) = 1;
			break;
		}
		case 12:
		{
			message.real()(i, 0) = 1;
			message.imag()(i, 0) = -3;
			break;
		}
		case 13:
		{
			message.real()(i, 0) = 1;
			message.imag()(i, 0) = -1;
			break;
		}
		case 14:
		{
			message.real()(i, 0) = 1;
			message.imag()(i, 0) = 3;
			break;
		}
		case 15:
		{
			message.real()(i, 0) = 1;
			message.imag()(i, 0) = 1;
			break;
		}
		}
	}
	return message;
}
Eigen::MatrixXcd& Message::DeMod16QAM()
{ 
	double temp = rows * cols;
	message = message.reshaped<double, double>(temp, 1);
	cols = message.cols();
	rows = message.rows();
	for (int i = 0; i < temp; i++)
	{
		/*0*/	if ((message.real()(i, 0) <= -2) && (message.imag()(i, 0) <= -2))
		{
			message.real()(i, 0) = 0;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*1*/if ((message.real()(i, 0) <= -2) && (message.imag()(i, 0) > -2) && (message.imag()(i, 0) <= 0))
		{
			message.real()(i, 0) = 1;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*2*/if ((message.real()(i, 0) <= -2) && (message.imag()(i, 0) > 2))
		{
			message.real()(i, 0) = 2;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*3*/if ((message.real()(i, 0) <= -2) && (message.imag()(i, 0) > 0) && (message.imag()(i, 0) <= 2))
		{
			message.real()(i, 0) = 3;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*4*/if ((message.real()(i, 0) > -2) && (message.real()(i, 0) <= 0) && (message.imag()(i, 0) <= -2))
		{
			message.real()(i, 0) = 4;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*5*/if ((message.real()(i, 0) > -2) && (message.real()(i, 0) <= 0) && (message.imag()(i, 0) > -2) && (message.imag()(i, 0) <= 0))
		{
			message.real()(i, 0) = 5;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*6*/if ((message.real()(i, 0) > -2) && (message.real()(i, 0) <= 0) && (message.imag()(i, 0) > 2))
		{
			message.real()(i, 0) = 6;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*7*/if ((message.real()(i, 0) > -2) && (message.real()(i, 0) <= 0) && (message.imag()(i, 0) > 0) && (message.imag()(i, 0) <= 2))
		{
			message.real()(i, 0) = 7;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*8*/if ((message.real()(i, 0) > 2) && (message.imag()(i, 0) < -2))
		{
			message.real()(i, 0) = 8;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*9*/if ((message.real()(i, 0) > 2) && (message.imag()(i, 0) > -2) && (message.imag()(i, 0) <= 0))
		{
			message.real()(i, 0) = 9;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*10*/if ((message.real()(i, 0) > 2) && (message.imag()(i, 0) > 2))
		{
			message.real()(i, 0) = 10;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*11*/if ((message.real()(i, 0) > 2) && (message.imag()(i, 0) > 0) && (message.imag()(i, 0) <= 2))
		{
			message.real()(i, 0) = 11;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*12*/if ((message.real()(i, 0) > 0) && (message.real()(i, 0) <= 2) && (message.imag()(i, 0) <= -2))
		{
			message.real()(i, 0) = 12;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*13*/if ((message.real()(i, 0) > 0) && (message.real()(i, 0) <= 2) && (message.imag()(i, 0) > -2) && (message.imag()(i, 0) <= 0))
		{
			message.real()(i, 0) = 13;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*14*/if ((message.real()(i, 0) > 0) && (message.real()(i, 0) <= 2) && (message.imag()(i, 0) > 2))
		{
			message.real()(i, 0) = 14;
			message.imag()(i, 0) = 0;
			continue;
		}
		/*15*/if ((message.real()(i, 0) > 0) && (message.real()(i, 0) <= 2) && (message.imag()(i, 0) > 0) && (message.imag()(i, 0) <= 2))
		{
			message.real()(i, 0) = 15;
			message.imag()(i, 0) = 0;
			continue;
		}

	}
	return message;
}
Eigen::MatrixXcd& Message::Reshape(int n_rows, int n_cols) // Аналогичен в Матлаб
{
	cols = n_cols;
	rows = n_rows;
	message = message.reshaped<double,double>(rows, cols);

	return message;
}
Eigen::MatrixXcd& Message::DFT_inv(std::string duration, int normalization_value) //обратное преобразование Фурье, с использованием временных векторов.
{    
	if (duration == "ROWS")
	{
		Eigen::VectorXcd temp(message.rows());
		for (int i = 0; i < cols; i++)
		{
			temp = message.col(i);
			plan_inv = fftw_plan_dft_1d(rows, (fftw_complex*)& temp(0), (fftw_complex*)& temp(0), FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(plan_inv);
			message.col(i) = temp/normalization_value;
			
		}
	}
	else
	{
		Eigen::VectorXcd temp(message.cols());
		for (int i = 0; i < cols; i++)
		{
			temp = message.row(i);
			plan_inv = fftw_plan_dft_1d(cols, (fftw_complex*)& temp(0), (fftw_complex*)& temp(0), FFTW_BACKWARD, FFTW_ESTIMATE);
			fftw_execute(plan_inv);
			message.row(i) = temp/normalization_value;
			
		}
	}
	
	return message;
}
Eigen::MatrixXcd& Message::DFT_fwd(std::string duration, int normalization_value) //Прямое преобразование.
{
	if (duration == "ROWS")
	{
		Eigen::VectorXcd temp(message.rows());
		for (int i = 0; i < cols; i++)
		{
			temp = message.col(i);
			plan_fwd = fftw_plan_dft_1d(rows, (fftw_complex*)& temp(0), (fftw_complex*)& temp(0), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(plan_fwd);
			message.col(i) = temp/normalization_value;
		
		}
	}
	else
	{
		Eigen::VectorXcd temp(message.cols());
		for (int i = 0; i < cols; i++)
		{
			temp = message.row(i);
			plan_fwd = fftw_plan_dft_1d(cols, (fftw_complex*)& temp(0), (fftw_complex*)& temp(0), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(plan_fwd);
			message.row(i) = temp/normalization_value;
		
		}
	}
	return message;
}
//Eigen::MatrixXcd& Message::AWGN(int SNR_dB)
//{
//	double power = 0;
//
//	for (int i = 0; i < rows; i++)
//	{
//		for (int j = 0; j < cols; j++)
//		{
//			power = power + (pow(message.real()(i, j), 2) + pow(message.imag()(i, j), 2));
//		}
//	}
//	power = message.size();
//	double snr = pow(10, SNR_dB / 10);
//	double noise_power = power / snr;
//	double sigma = sqrt(noise_power);
//	Eigen::MatrixXcd temp_noise(rows, cols);
//	double** rand = randn(rows, cols);
//	for (int i = 0; i < rows; i++)
//	{
//		for (int j = 0; j < cols; j++)
//		{
//			temp_noise.real()(i, j) = sqrt(noise_power / 2) * rand[i][j];
//			temp_noise.imag()(i, j) = rand[i][j];
//		}
//
//	}
//
//	for (int i = 0; i < rows; i++)
//	{
//		for (int j = 0; j < cols; j++)
//		{
//			message.real()(i, j) = message.real()(i, j) + temp_noise.real()(i, j);
//			message.imag()(i, j) = message.imag()(i, j) + temp_noise.imag()(i, j);
//		}
//	}
//
//	for (int i = 0; i < rows; i++)
//	{
//		delete[] rand[i];
//	}
//	delete[] rand;
//
//
//	return message;
//}
void Message::ShowMessage() // Вывод матрицы сообщения
{
	std::cout << "Message now: \n" << message << "\n\n";
};
int Message::Size() //Взятие размера сообщения
{
	return message.size();
};
int Message::Rows() // Взятие кол-ва строк
{
	return rows;
}
int Message::Cols() //Взятие столбцов
{
	return cols;
}
double Message::GetReal(int i, int j)  // Взятие реальной части
{
	return message.real()(i, j);
}
double Message::GetImag(int i, int j) // Взятие мнимой части.
{
	return message.imag()(i, j);
}
Eigen::MatrixXcd& Message::DeModBPSK_EVK()
{
	int size = rows * cols;
	message = message.reshaped<int, int>(size, 1);
	Eigen::MatrixXd::Index minRow, minCol;
	Eigen::MatrixXcd number(2, 2);
	number(0, 0) = std::complex<double>(-1.0, 0.0);
	number(0, 1) = std::complex<double>(1.0, 0.0);
	number(1, 0) = std::complex<double>(1.0, 0.0);
	number(1, 1) = std::complex<double>(0.0, 0.0);
	Eigen::VectorXd temp(2);
	
	for (int i = 0; i < size; i++)
	{
		for (int m = 0; m < 2; m++)
		{	
			 temp(m) = sqrt(pow(message.real()(i, 0) - number.real()(0, m), 2) +  pow(message.imag()(i, 0) - number.imag()(0, m), 2));
		}
		temp.minCoeff(&minRow, &minCol);
		message(i) = number(1, minRow);
	}
	cols = message.cols();
	rows = message.rows();
	return message;
};
Eigen::MatrixXcd& Message::DeModQPSK_EVK() 
{
	int size = rows * cols;
	double value = 0.70710548251124;
	message = message.reshaped<int, int>(size, 1);
	Eigen::MatrixXd::Index minRow, minCol;
	Eigen::MatrixXcd number(2, 4);
	number(0, 0) = std::complex<double>(value, value);
	number(0, 1) = std::complex<double>(value, -1 * value);
	number(0, 2) = std::complex<double>(-1 * value, value);
	number(0, 3) = std::complex<double>(-1 * value, -1 * value);
	number(1, 0) = std::complex<double>(0.0, 0.0);
	number(1, 1) = std::complex<double>(1.0, 0.0);
	number(1, 2) = std::complex<double>(2.0, 0.0);
	number(1, 3) = std::complex<double>(3.0, 0.0);
	
	Eigen::VectorXd temp(4);
	for (int i = 0; i < size; i++)
	{
		for (int m = 0; m < 4; m++)
		{
			temp(m) = sqrt(pow(message.real()(i, 0) - number.real()(0, m), 2) + pow(message.imag()(i, 0) - number.imag()(0, m), 2));
		}
		temp.minCoeff(&minRow, &minCol);
		message(i, 0) = number(1, minRow);
	}
	cols = message.cols();
	rows = message.rows();
	return message;

};
Eigen::MatrixXcd& Message::DeMod16QAM_EVK()
{
	int size = rows * cols;
	message = message.reshaped<int, int>(size, 1);
	Eigen::MatrixXd::Index minRow, minCol;
	std::complex<double> values[16]{ {-3.0, -3.0}, {-3.0, -1.0}, {-3.0, 3.0}, {-3.0,1.0}, {-1.0,-3.0}, {-1.0,-1.0},
	{-1.0, 3.0}, {-1.0, 1.0}, {3.0, -3.0}, {3.0, -1.0}, {3.0, 3.0}, {3.0, 1.0}, {1.0, -3.0}, {1.0, -1.0}, {1.0, 3.0},
	{1.0, 1.0} };
	Eigen::MatrixXcd number(2, 16);
	int i = 0; double j = 0.0;
	for (; i < 16 ; i++, j++)
	{   
		number(0, i) = values[i];
		number(1, i) = std::complex<double>(j, 0.0);
	}
	Eigen::VectorXd temp(16);
	for (int i = 0; i < size; i++)
	{
		for (int m = 0; m < 16; m++)
		{
			temp(m) = sqrt(pow(message.real()(i, 0) - number.real()(0, m), 2) + pow(message.imag()(i, 0) - number.imag()(0, m), 2));
		}
		temp.minCoeff(&minRow, &minCol);
		message(i, 0) = number(1, minRow);
	}
	cols = message.cols();
	rows = message.rows();
	return message;

};
Eigen::MatrixXcd& Message::IFFTSHIFT()
{
	Eigen::MatrixXcd temp(message.rows(), message.cols());
	temp = circShift(message, ceil(message.rows() / 2), ceil(message.cols() / 2));
	message = temp;
	return message;
}
Eigen::MatrixXcd& Message::IFFTSHIFT(int dir)
{
	Eigen::MatrixXcd temp(message.rows(), message.cols());
		if (dir == 1)
		{
			temp = circShift(message, ceil(message.rows()/2), 0);
			message = temp;
		}
		else 
		{
			temp = circShift(message, ceil(message.cols() / 2), 0);
			message = temp;
		}
	return message;
};
Eigen::MatrixXcd& Message::FFTSHIFT()
{
	Eigen::MatrixXcd temp(message.rows(), message.cols());
	temp = circShift(message, floor(message.rows() / 2), floor(message.cols() / 2));
	message = temp;
	return message;
};
Eigen::MatrixXcd& Message::FFTSHIFT(int dir)
{
	Eigen::MatrixXcd temp(message.rows(), message.cols());
	if (dir == 1)
	{
		temp = circShift(message, floor(message.rows() / 2), 0);
		message = temp;
	}
	else
	{
		temp = circShift(message, 0, floor(message.cols() / 2));
		message = temp;
	}
	return message;
};
void Message::getFeatures_14dwt3stat5cum(Eigen::MatrixXd& features)
	{
		int p = 1;
		for (int i = 0; i < 22; i++)
		{
			features(i) = 0;
		}
		Eigen::MatrixXd cA_real; Eigen::MatrixXd cA_imag;
		Eigen::MatrixXd cD_real; Eigen::MatrixXd cD_imag;
		Eigen::MatrixXd Real_Signal(message.rows(), message.cols());
		Eigen::MatrixXd Imag_Signal(message.rows(), message.cols());
		for (int i = 0; i < cols; i++)
		{
			Real_Signal.col(i) = message.col(i).real();
			Imag_Signal.col(i) = message.col(i).imag();
		}
		Real_Signal = Real_Signal.reshaped(Real_Signal.rows() * Real_Signal.cols(), 1);
		Imag_Signal = Imag_Signal.reshaped(Imag_Signal.rows() * Imag_Signal.cols(), 1);
		std::string z[7]{ "db1","db2","db3","db4", "db5", "db6", "db7" };
		int* ptrr;
		for (int m = 0, j = 7; m < 7; m++, j++)
		{
			const char* ptr = z[m].c_str();
			DWT(Real_Signal, ptr, cA_real, cD_real);
			DWT(Imag_Signal, ptr, cA_imag, cD_imag);
			int size = cA_real.rows();
			ptrr = &size;
			Eigen::MatrixXcd cA(size, 1); Eigen::MatrixXcd cD(size, 1);
			for (int i = 0; i < size; i++)
			{
				cA(i, 0) = { cA_real(i),cA_imag(i) };
				cD(i, 0) = { cD_real(i),cD_imag(i) };
			}
			for (int i = 0; i < size; i++)
			{
				features(m) = features(m) + pow(sqrt((pow(cA(i,0).real(), 2) + pow(cA(i,0).imag(), 2))), p);	
				features(j) = features(j) + pow(sqrt((pow(cD(i,0).real(), 2) + pow(cD(i,0).imag(), 2))), p);
			}
			features(m) = features(m) / *ptrr;
			features(j) = features(j) / *ptrr;
		}

		int N_FFT = pow(2, 14);
		int size = Real_Signal.rows();
		Eigen::MatrixXcd Hilbert_signal(size, 1);
		for (int i = 0; i < size; i++)
		{
			Hilbert_signal(i, 0) = {Real_Signal(i,0), 0};
		}
		Hilbert(Hilbert_signal, Hilbert_signal);
		for (int i = 0; i < size; i++)
		{
			Hilbert_signal(i) = { sqrt((pow(Hilbert_signal(i).real(), 2) + pow(Hilbert_signal(i).imag(), 2))), 0 };
		}
		double m_a = 0;
		for (int i = 0; i < size; i++)
		{
			m_a += Hilbert_signal(i).real();
		}
		m_a /= size;
		Hilbert_signal /= m_a;
		std::complex<double> value(1, 0);
		for (int i = 0; i < size; i++)
		{
			Hilbert_signal(i) = Hilbert_signal(i) - value;
		}
		
		Eigen::MatrixXcd value_arr(size,1);
		fftw_plan plan = fftw_plan_dft_1d(size, (fftw_complex*)& Hilbert_signal(0), (fftw_complex*)& value_arr(0), FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		value_arr = circShift(value_arr, floor(size / 2), 0);
		for (int i = 0; i < size; i++)
		{
			value_arr(i) = { sqrt((pow(value_arr(i).real(), 2) + pow(value_arr(i).imag(), 2))), 0 };
		}
		double gamma_max = pow(value_arr.real().maxCoeff(), 2);
		gamma_max /= size;
		double L = 0, value_L = 0;
		for (int i = 0; i < size; i++)
		{
			L += pow(Hilbert_signal(i).real(), 2);
			value_L += abs(Hilbert_signal(i).real());
		}
		L /= size; value_L /= size; value_L = pow(value_L, 2);
		double sigma_aa = sqrt(abs(L - value_L));
		L = log10(abs(L - value_L));

		std::complex<double> c20 = { 0,0 }, c40 = { 0,0 }, c41 = { 0,0 }; double c21 = 0, c42 = 0;
		std::complex<double> value_c = { 0,0 };
		value_arr.resize(rows, cols);
		value_arr = message;
		value_arr = value_arr.reshaped(size, 1);

		for (int i = 0; i < size; i++)
		{
			value_c = { pow(value_arr(i).real(),2)- pow(value_arr(i).imag(),2), 2 * value_arr(i).imag() * value_arr(i).real() };
			c20 = value_c + c20;
		}
		c20 = { c20.real() / size,c20.imag() / size };

		for (int i = 0; i < size; i++)
		{
			c21 = c21 +(pow(value_arr(i).real(), 2) + pow(value_arr(i).imag(), 2));
		
		}
		c21 /= size;

		for (int i = 0; i < size; i++)
		{
			value_c = {pow(value_arr(i).real(),4) + pow(value_arr(i).imag(),4) - 6 * pow(value_arr(i).real() * value_arr(i).imag(),2),
			4 * (value_arr(i).imag() * pow(value_arr(i).real(),3) - value_arr(i).real() * pow(value_arr(i).imag(),3))};
			c40 = value_c + c40;
		}
		c40 = { c40.real() / size,c40.imag() / size }; c40 = {c40.real() - 3 * c20.real(),c40.imag() - 3 * c20.imag()};

		for (int i = 0; i < size; i++)
		{
			value_c = { pow(value_arr(i).real(), 3) - 3 * value_arr(i).real() * pow(value_arr(i).imag(),2),
			3 * pow(value_arr(i).real(), 2) * value_arr(i).imag() -  pow(value_arr(i).imag(), 3)};
			value_c = { value_c.real() * value_arr(i).real() + value_c.imag() * value_arr(i).imag() };
			c41 = value_c + c41;
		}
		c41 = { c41.real() / size,c41.imag() / size }; 
		c41 = { c41.real() - 3 * c20.real() * c21 ,c41.imag() - 3 * c20.imag() * c21};

		for (int i = 0; i < size; i++)
		{
			c42 += pow((pow(value_arr(i).real(),2) + pow(value_arr(i).imag(), 2)),2);
		}
		c42 = c42 / size - (pow(c20.real(), 2) + pow(c20.imag(), 2)) - 2 * pow(c21, 2);

		features(14) = abs(gamma_max);
		features(15) = abs(L);
		features(16) = abs(sigma_aa);
		features(17) = sqrt(pow(c20.real(), 2) + pow(c20.imag(), 2));
		features(18) = abs(c21);
		features(19) = sqrt(pow(c40.real(), 2) + pow(c40.imag(), 2));
		features(20) = sqrt(pow(c41.real(), 2) + pow(c41.imag(), 2));
		features(21) = abs(c42);
	};

