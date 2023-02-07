#include "Purelin.h"

Eigen::MatrixXd& purelin(Eigen::MatrixXd& Weigth, Eigen::MatrixXd& features, Eigen::MatrixXd& b)
{
	int col = features.cols(), row = Weigth.rows();

	Eigen::MatrixXd temp(row, col);
	temp = Weigth * features;
	for (int i = 0; i < col; i++)
	{
		temp.col(i) = temp.col(i) + b;
	}
	return temp;
}
