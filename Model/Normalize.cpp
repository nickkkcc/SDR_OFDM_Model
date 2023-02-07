#include "Normalize.h"

void normalize(Eigen::MatrixXd& Input_features, Eigen::MatrixXd& Output_features, Net_Info& Info)
{
	for (int i = 0; i < Input_features.cols(); i++)
	{
			Input_features.col(i) = Input_features.col(i) -  Info.offset;
	}
	Input_features *= Info.gain;
	for (int i = 0; i < Input_features.cols(); i++)
	{
		for (int j = 0; j < Input_features.rows(); j++)
		{
			Input_features.col(i)(j) = Input_features.col(i)(j) * Info.gain(j);
		}
	}
	for (int i = 0; i < Input_features.size(); i++)
	{
		for (int j = 0; j < Input_features.cols(); j++)
		{
			Input_features(i, j) += Info.ymin;
		}
	}
	Output_features = Input_features;
}