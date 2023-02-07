#include<iostream>
#include <math.h>


double RandomDouble(int min, int max, int precision) //Функция создания равномерных значений в промежутке с определенной точностью.
{
	double value = rand() % (int)pow(10, precision);
	return(min + (value / pow(10, precision)) * (max - min));

};

double** randn(int rows, int cols) //Создание массива стандартно - отклоненных величин с помощью метода Бокса-Мюллера.
{

	double** arr = new double* [rows];

	for (int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			double u = 0;
			double v = 0;
			while (((pow(u, 2) + pow(v, 2)) >= 1) || ((pow(u, 2) + pow(v, 2)) <= 0))
			{
				u = RandomDouble(-2, 2, 3);
				v = RandomDouble(-2, 2, 3);
			}
			double s = pow(u, 2) + pow(v, 2);
			arr[i][j] = u * sqrt((-2 * log(s)) / s);

		}
	}
	return arr;
}