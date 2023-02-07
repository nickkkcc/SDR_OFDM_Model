#include "Run_Hardware.h"
#include "tansig.h"
#include "Purelin.h"
void Run_Hardware(Eigen::MatrixXd& features, Eigen::MatrixXd& Y, Net_Info& Info) 
{
	int sample = features.cols();
    normalize(features, features, Info);
	Y.resize(3, sample);
	Eigen::MatrixXd Y1(Info.IW.rows(), sample);
	Eigen::MatrixXd Y2(Info.LW2.rows(), sample); 
	Eigen::MatrixXd Y3(Info.LW3.rows(), sample);
	Eigen::MatrixXd Y4(Info.LW4.rows(), sample);
	
	Y1 = tansig(Info.IW, features, Info.b1);
	Y2 = tansig(Info.LW2, Y1, Info.b2);
	Y3 = tansig(Info.LW3, Y2, Info.b3);
	Y4 = purelin(Info.LW4, Y3, Info.b4);
	Y = Y4;
};
