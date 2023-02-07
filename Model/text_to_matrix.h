#pragma once
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <string>

void getTextToMatrix(std::string&, Eigen::MatrixXd&);
void getTextToMatrix(std::string* PATH, double&);
void getTextToMatrix(std::string& PATH_real, std::string& PATH_imag, Eigen::MatrixXcd&);
