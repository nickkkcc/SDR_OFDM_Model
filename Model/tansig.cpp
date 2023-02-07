#include "tansig.h"

Eigen::MatrixXd& tansig(Eigen::MatrixXd& Weigth, Eigen::MatrixXd& features, Eigen::MatrixXd& b)
{
	int col = features.cols(), row = Weigth.rows();

	Eigen::MatrixXd temp(row, col);
	temp = Weigth * features;
	for (int i = 0; i < col; i++)
	{
		temp.col(i) = temp.col(i) + b;
	}
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			temp(i,j) = (2 / (1 + exp(-2 * temp(i,j))) - 1);
		}
	}
	return temp;

};
