#include "SourceVolume.h"
#include <vector>
#include <string>
#include <iostream>
#include <valarray>
#include <fstream>

double getExpectation(vector<int>& x)
{
	double sum = 0;
	for (auto p : x)
	{
		sum += p;
	}
	sum /= x.size();
	return sum;
}

double suma(vector<double> a)
{
	double s = 0;
	for (int i = 0; i < a.size(); i++)
	{
		s = s + a[i];
	}
	return s;
}

double mean(vector<double>& a)
{
	return suma(a) / a.size();
}

double sqsum(vector<double>& a)
{
	double s = 0;
	for (int i = 0; i < a.size(); i++)
	{
		s = s + pow(a[i], 2);
	}
	return s;
}

double stdev(vector<double>& nums)
{
	double N = nums.size();
	return pow(sqsum(nums) / N - pow(suma(nums) / N, 2), 0.5);
}

vector<double> operator-(vector<double>& a, double b)
{
	vector<double> retvect;
	for (int i = 0; i < a.size(); i++)
	{
		retvect.push_back(a[i] - b);
	}
	return retvect;
}

vector<double> operator*(vector<double> a, vector<double> b)
{
	vector<double> retvect;
	for (int i = 0; i < a.size(); i++)
	{
		retvect.push_back(a[i] * b[i]);
	}
	return retvect;
}

double pearsoncoeff(vector<double> X, vector<double> Y)
{
	return suma((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}
double getLLC(const vector<vector<unsigned char>>& regular_data, int index, const int dimension_x, const int dimension_y, const int dimension_z, int v_idx_a, int v_idx_b, const int window_size)
{
	const auto k = index / (dimension_x*dimension_y);
	const auto i = index % dimension_x;
	const auto j = (index % (dimension_x*dimension_y)) / dimension_x;

	//针对每个体素位置，获取窗口内的局部相关性 EX, EY, EXY, EXEY, EX2, (EX)2, EY2, (EY)2

	vector<double> X, X2, Y, Y2, XY;
	for (auto zvalue = k - window_size; zvalue <= k + window_size; zvalue++) {
		if (zvalue < 0 || zvalue >= dimension_z) continue;
		//for y
		for (auto yvalue = j - window_size; yvalue <= j + window_size; yvalue++) {
			if (yvalue < 0 || yvalue >= dimension_y) continue;
			//for x
			for (auto xvalue = i - window_size; xvalue <= i + window_size; xvalue++) {
				if (xvalue < 0 || xvalue >= dimension_x) continue;
				const auto x_buf = regular_data[v_idx_a][zvalue*dimension_x*dimension_y + yvalue * dimension_x + xvalue];
				const auto y_buf = regular_data[v_idx_b][zvalue*dimension_x*dimension_y + yvalue * dimension_x + xvalue];
				X.push_back(x_buf);
				//X2.push_back(x_buf*x_buf);
				Y.push_back(y_buf);
				//Y2.push_back(y_buf*y_buf);
				//XY.push_back(x_buf*y_buf);
			}
		}
	}

	auto p = pearsoncoeff(X, Y);
	//cout << p << endl;
	if (-1 < p &&p < 1) return abs(p);
	else return 0;

	auto EX = 0.0, EY = 0.0, EX2 = 0.0, EY2 = 0.0, EXY = 0.0;;
	for (auto a = 0; a < X.size(); a++)
	{
		EX += X[a];
		EY += Y[a];
		EX2 += X[a] * X[a];
		EY2 += Y[a] * Y[a];
		EXY += X[a] * Y[a];
	}
	for (auto a = 0; a < X.size(); a++)
	{
		EX /= X.size();
		EY /= X.size();
		EX2 /= X.size();
		EY2 /= X.size();
		EXY /= X.size();
	}

	auto cov_xy = EXY - EX * EY;
	auto DX = sqrt(EX2 - EX * EX);
	auto DY = sqrt(EY2 - EY * EY);

	//两者方差均为0
	if (DX < FLT_EPSILON&&DY < FLT_EPSILON)
		return 1;

	if (DX < FLT_EPSILON) DX = 1;
	if (DY < FLT_EPSILON) DY = 1;

	//if(abs(cov_xy / sta_x * sta_y)>1000)
	{
		//int a = 12000;
		//cout << cov_xy / DX * DY <<"\t"<< cov_xy << "\t" << EX << "\t" << EY << "\t" << EX2 << "\t" << EY2 << endl;
		//getchar();
	}

	return min(abs(cov_xy / DX * DY), 1.0);
}


double getLLC2(const vector<vector<unsigned char>>& regular_data, int index, const int dimension_x, const int dimension_y, const int dimension_z, int v_idx_a, int v_idx_b, const int window_size)
{
	const auto k = index / (dimension_x*dimension_y);
	const auto i = index % dimension_x;
	const auto j = (index % (dimension_x*dimension_y)) / dimension_x;

	//针对每个体素位置，获取窗口内的局部相关性 EX, EY, EXY, EXEY, EX2, (EX)2, EY2, (EY)2

	vector<double> X, X2, Y, Y2, XY;

	float hist1D1[256];
	float hist1D2[256];
	float hist2D[256 * 256];
	//Init histogram
	for (int i = 0; i < 256; i++) {
		hist1D1[i] = 0;
		hist1D2[i] = 0;
	}
	for (int i = 0; i < 256 * 256; i++) {
		hist2D[i] = 0;
	}
	auto cnt = 0;
	for (auto zvalue = k - window_size; zvalue <= k + window_size; zvalue++) {
		if (zvalue < 0 || zvalue >= dimension_z) continue;
		//for y
		for (auto yvalue = j - window_size; yvalue <= j + window_size; yvalue++) {
			if (yvalue < 0 || yvalue >= dimension_y) continue;
			//for x
			for (auto xvalue = i - window_size; xvalue <= i + window_size; xvalue++) {
				if (xvalue < 0 || xvalue >= dimension_x) continue;
				const auto x_buf = regular_data[v_idx_a][zvalue*dimension_x*dimension_y + yvalue * dimension_x + xvalue];
				const auto y_buf = regular_data[v_idx_b][zvalue*dimension_x*dimension_y + yvalue * dimension_x + xvalue];
				X.push_back(x_buf);
				//X2.push_back(x_buf*x_buf);
				Y.push_back(y_buf);
				//Y2.push_back(y_buf*y_buf);
				//XY.push_back(x_buf*y_buf);

				hist1D1[x_buf]++;
				hist1D2[y_buf]++;
				hist2D[x_buf * 256 + y_buf]++;
				cnt++;
			}
		}
	}
	

	double EX = 0.0f, EY = 0.0f, EXY = 0.0f;
	double DX = 0.0f, Ex2 = 0.0f;
	double DY = 0.0f, Ey2 = 0.0f;

	//Calculate EX, EY, Ex2, Ey2
	for (int i = 0; i < 256; i++) {
		hist1D1[i] /= X.size();
		Ex2 += i * i*hist1D1[i];
		EX += i * hist1D1[i];

		hist1D2[i] /= X.size();
		Ey2 += i * i*hist1D2[i];
		EY += i * hist1D2[i];
	}
	DX = Ex2 - EX * EX;
	DY = Ey2 - EY * EY;

	if (abs(EXY - EX * EY) < FLT_EPSILON)
		return 0;

	//两者方差均为0
	if (DX < FLT_EPSILON&&DY < FLT_EPSILON)
		return 1;

	if (DX < FLT_EPSILON) DX = 1;
	if (DY < FLT_EPSILON) DY = 1;

	//Calculate EXY
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			hist2D[i * 256 + j] /= X.size();
			EXY += i * j*hist2D[i * 256 + j];
		}
	}
	//Get person correlation coefficient and detect the error of out of range
	double pearson = (EXY - EX * EY) / (sqrt(DX)*sqrt(DY));
	//if (pearson > 1 - FLT_EPSILON || pearson < -1 + FLT_EPSILON) {
	//	cout << "Result Error: Pearson correlation coefficient!" << endl;
	//}
	//cout << pearson << "\t" << EXY << endl;
	//Need return absolute value?
	return abs(pearson);

}


vector<double> calcLocalCorrelationCoefficient(const vector<vector<unsigned char>>& regular_data, const vector<int>& index_array, 
	const int dimension_x, const int dimension_y, const int dimension_z, const int window_size)
{
	vector<double> llc_array(dimension_x*dimension_y*dimension_z,1);
	if(index_array.size() >= 2)
	{
		for(auto i=0;i<index_array.size();i++)
		{
			for(auto j=i+1;j<index_array.size();j++)
			{
#pragma omp parallel for
				for (auto index = 0; index < regular_data[i].size(); index++)
				{
					llc_array[index] = min(llc_array[index],
						getLLC(regular_data, index, dimension_x, dimension_y, dimension_z, index_array[i], index_array[j], window_size));
				}
			}
		}
	}
	cout << "Local correlation coefficient has been calculated." << endl;
	return llc_array;
}

void saveVolume(vector<unsigned char>& cs, const vector<double>& xes, double similarity_threshold)
{
	ofstream mask_file("./result/mask.raw", ios_base::binary);
	unsigned char buf = 0;
	for (auto i = 0; i < cs.size(); i++)
	{
		if(xes[i]>=similarity_threshold)
			mask_file.write((char*)&cs[i], sizeof(cs[i]));
		else
			mask_file.write((char*)&buf, sizeof(cs[i]));
	}
	mask_file.close();
	std::cout << "The correlation file has been saved." << std::endl;
}

int main(){

	vector<string> file_list;
	file_list.push_back("F:\\atmosphere\\timestep21_float\\Pf21.bin");
	file_list.push_back("F:\\atmosphere\\timestep21_float\\TCf21.bin");
	file_list.push_back("F:\\atmosphere\\timestep21_float\\QVAPORf21.bin");
	file_list.push_back("F:\\atmosphere\\timestep21_float\\CLOUDf21.bin");
	file_list.push_back("F:\\atmosphere\\timestep21_float\\PRECIPf21.bin");

	file_list.push_back("F:\\atmosphere\\timestep21_float\\QVAPORf16.bin");
	file_list.push_back("F:\\atmosphere\\timestep21_float\\TCf16.bin");

	int dimension_x = 500, dimension_y = 500, dimension_z = 100;

	SourceVolume volume(file_list, dimension_x, dimension_y, dimension_z, "float", 256, 256);
	
	volume.loadVolume();
	volume.loadRegularVolume();

	vector<vector<unsigned char>> regular_data(file_list.size());
	for(auto i=0;i<regular_data.size();i++)
	{
		regular_data[i] = *(volume.getRegularVolume(i));
	}

	auto llc_array = calcLocalCorrelationCoefficient(regular_data, {2, 3, 4}, dimension_x, dimension_y, dimension_z, 3);

	auto similarity_threshold = 0.5;
	//vector<double> llc_array;
	saveVolume(regular_data[1], llc_array, similarity_threshold);

	getchar();
}
