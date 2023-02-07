#include <iostream>
#include "ErrorCount.h"
#include "MessageClass.h"


double ErrorCount(Message& A, Message& B) // Вычисление битовой ошибки
{
	double error = 0;

	if (A.Rows() != B.Rows() || A.Cols() != B.Cols())
	{
		std::cout << "Dimensions of messages are not equal!\n";
		return -1;
	}
	else
	{
		for (size_t i = 0; i < A.Rows(); i++)
		{
			for (int j = 0; j < A.Cols(); j++)
			{
				if (A.GetReal(i, j) != B.GetReal(i, j))
				{
					error++;
				}
			}
		}
	}
	return error;
};