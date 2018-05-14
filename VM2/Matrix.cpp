#pragma once
#include<iostream>
#include<math.h>
#include<vector>
#include<string>
#include"Matrix.h"
#include"func.h"
Matrix::Matrix(int M, int N):dimM(M), dimN(N)
{
	mas = new double*[dimM];
	for (int i = 0; i < dimM; i++)
		mas[i] = new double[dimN];
	for (int i = 0; i < dimM; i++)
		for (int j = 0; j < dimN; j++)
			mas[i][j] = 0;
}
//Matrix::~Matrix()
//{
//delete[] mas;
//}

void Matrix::show()const
{

	for (int i = 0; i < dimM; i++)
	{
		std::cout << std::endl;
		for (int j = 0; j < dimN; j++)
		{
			if (abs(mas[i][j])<E)
				std::cout << 0 << " ";
			else
				std::cout << mas[i][j] << " ";
		}
	}
	std::cout << std::endl << std::endl;
}
bool Matrix::operator==(const Matrix &A)
{


	for (int i = 0; i < dimM; i++)
		for (int j = 0; j < dimN; j++)
			if (abs(mas[i][j] - A.mas[i][j]) > E)
			{

				return false;
			}
	return true;

}
void Matrix::randomDP()
{
	for (int i = 0; i < dimM; i++)
	{
		double sum = 0;
		for (int j = 0; j < dimN; j++)
		{
			mas[i][j] = rand() % 100 - 50;
			if(i!=j)
				sum += abs(mas[i][j]);
		}
		mas[i][i] =abs(mas[i][i])+ sum;
	}
}
const Matrix Matrix::operator=(const Matrix &A)
{
	for (int i = 0; i < A.dimM; i++)
		for (int j = 0; j < A.dimN; j++)
			mas[i][j] = A.mas[i][j];
	return *this;
}
const Matrix Matrix::operator=(const double inMas[])
{
	for (int i=0;i<dimM;++i)
		for (int j = 0; j < dimN; ++j)
		{
			mas[i][j] = inMas[i*dimN + j];
		}
	
	return *this;
}
double* Matrix::operator[](int i)const
{
	return mas[i];
}
Matrix Matrix::operator*(const Matrix &A)
{
	Matrix temp(dimM, A.dimN);
	for (int i = 0; i < dimM; i++)
		for (int j = 0; j < A.dimN; j++)
			for (int k = 0; k < A.dimM; k++)
			{


				temp[i][j] += mas[i][k] * A[k][j];
			}
	return temp;
}
Matrix Matrix::operator*(double K)
{
	Matrix temp(dimM, dimN);
	for (int i = 0; i < dimM; i++)
		for (int j = 0; j < dimN; j++)
			temp[i][j] =mas[i][j]* K;
	return temp;
}
Matrix Matrix::operator+(const Matrix & A)
{
	Matrix temp(dimM, dimN);
	for (int i = 0; i < dimM; i++)
		for (int j = 0; j < dimN; j++)
			temp[i][j] = mas[i][j] +A[i][j];
	return temp;
}
//check!
Matrix Matrix::operator-(Matrix & A)
{
	return *this+ A*(-1.0);
}
void Matrix::SwapRows(int f, int s)
{
	if (f == s || f == -1)
		return;

	for (int i = 0; i < dimN; i++)
	{
		std::swap(mas[s][i], mas[f][i]);
	}
	sign = -sign;
	//std::cout << "swapR " << f << " " << s << "\n" << "sign " << sign << "\n";
}
void Matrix::SwapCol(int f, int s)
{
	if (f == s || f == -1)
		return;

	for (int i = 0; i < dimM; i++)
	{
		std::swap(mas[i][s], mas[i][f]);
	}
	sign = -sign;
	//std::cout << "swapC " << f << " " << s << "\n" << "sign " << sign << "\n";
}
Matrix Matrix::Transpose()
{
	Matrix tmp(dimM, dimN);
	tmp = *this;
	for (int i = 0; i < dimM; i++)
		for (int j = 0; j < i; j++)
		{
			std::swap(tmp[i][j], tmp[j][i]);
		}
	return tmp;
}
double Matrix::CondNum()
{
	return this->Norm()*inverse(*this).Norm();
}
double Matrix::Determinant()
{
	double det = sign;
	for (int i = 0; i < dimM; i++)
	{
		det *= mas[i][i];
	}

	return det;
}
int Matrix::Rank(Matrix &U)
{
	int rank = 0;
	for (int i = 0; i < U.dimM; i++)
	{
		if (U[i][i] != 0)
			rank++;
	}
	return rank;
}
double Matrix::Norm()
{
	double max = 0;
	for (int i = 0; i < dimM; i++)
	{
		double sum = 0;
		for (int j = 0; j < dimM; j++)
		{
			sum += abs(mas[i][j]);
		}
		if (max < sum)
			max = sum;
	}
	return max;
}