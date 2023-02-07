#include "text_to_matrix.h"

void getTextToMatrix(std::string& PATH, Eigen::MatrixXd& Output)
{
	std::ifstream file(PATH,std::ios::in);
	if (file.is_open())
	{
		std::cout << PATH << " is open!" << std::endl;
		int rows = -1, cols = 0;
		std::string str;
		while (!file.eof())
		{
			getline(file, str);
			rows++;
		}
		file.close();
		file.open(PATH);
		std::vector<double> values(rows);
		for (int i = 0; i < rows; i++)
		{
			file >> values[i];
		}
		/*std::copy(std::istream_iterator<double>(file), std::istream_iterator<double>(), std::back_inserter(values));*/
		cols = values.size() / rows;
		Output.resize(rows, cols);
		int xx, index = 0;
		for (int i = 0; i < rows; i++)
		{
			xx = 0;
			while (xx < cols)
			{
				Output.row(i)(xx) = values[index];
				xx++;
				index++;
			}
		}
		file.close();
	}
	else
	{
		std::cout << PATH << " isn't open!" << std::endl;
	}
}
void getTextToMatrix(std::string* PATH, double& Output)
{
	std::ifstream file(*PATH, std::ios::in);
	if (file.is_open())
	{
		std::cout << *PATH << " is open!" << std::endl;
		std::vector<double> values;
		std::copy(std::istream_iterator<double>(file), std::istream_iterator<double>(), std::back_inserter(values));
		file.close();
		Output = values[0];
	}
	else 
	{
		std::cout << PATH << " isn't open!" << std::endl;
	}

}
void getTextToMatrix(std::string& PATH_real, std::string& PATH_imag, Eigen::MatrixXcd& Output)
{
	std::ifstream file(PATH_real, std::ios::in);
	
	if (file.is_open())
	{
		std::cout << " is open!" << std::endl;
		int rows = -1, cols = 0;
		
		std::string str;
		while (!file.eof())
		{
			getline(file, str);
			rows++;
		}
		file.close();

		file.open(PATH_real);
		std::vector<double> values(rows);
		for (int i = 0; i < rows; i++)
		{
			file >> values[i];
		}
		/*std::copy(std::istream_iterator<double>(file), std::istream_iterator<double>(), std::back_inserter(values));*/
		cols = values.size() / rows;
		int rows1 = rows, cols1 = cols;
		Output.resize(32, rows/32);
		rows = Output.rows(); cols = Output.cols();
		int xx, index = 0;
		for (int i = 0; i < cols; i++)
		{
			xx = 0;
			while (xx < rows)
			{
				Output.real().col(i)(xx) = values[index];
				xx++;
				index++;
			}
		}
		file.close();
		std::ifstream file1(PATH_imag, std::ios::in);
		std::vector<double> values1(rows1);
		for (int i = 0; i < rows1; i++)
		{
			file1 >> values1[i];
		}
		/*std::copy(std::istream_iterator<double>(file1), std::istream_iterator<double>(), std::back_inserter(values1));*/
		int cc, indexx = 0;
		for (int i = 0; i < cols; i++)
		{
			cc = 0;
			while (cc < rows)
			{
				Output.imag().col(i)(cc) = values1[indexx];
				cc++;
				indexx++;
			}
		}
		file1.close();
	}
	else
	{
		std::cout << " isn't open!" << std::endl;
	}
}
