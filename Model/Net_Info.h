#pragma once
#include <Eigen/Dense>
class Net_Info
{
public:
	Eigen::MatrixXd gain;
	Eigen::MatrixXd offset;
	double ymin;
	Eigen::MatrixXd IW;
	Eigen::MatrixXd LW2;
	Eigen::MatrixXd LW3;
	Eigen::MatrixXd LW4;
	Eigen::MatrixXd b1;
	Eigen::MatrixXd b2;
	Eigen::MatrixXd b3;
	Eigen::MatrixXd b4;
	Eigen::MatrixXd signal;
	Eigen::MatrixXd features;

	Net_Info();


	~Net_Info();
};

