#pragma once 
#include<iostream>
#include"Matrix.h"
#include"func.h"
#include<functional>
#include<vector>
#include <time.h> 
const double Eps = 0.00000001;

auto Matrix_x02 = []()
{
	Matrix result(10, 1);
	double x[] =
	{
		0.5, 0.5, 1.5, -1.0, -0.2, 1.5, 0.5, -0.5, 1.5, -1.5
		//0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
	};
	return result = x;

};

auto Matrix_x0 = []()
{
	Matrix result(10, 1);
	double x[] =
	{
		0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5
		//0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
	};
	return result = x;

};
auto SNAE_matrix = [](const Matrix &x)
{
	Matrix result(10, 1);
	double SNAE_fm[] = {
		cos(x[0][0] * x[1][0]) - exp(-3. * x[2][0]) + x[3][0] * x[4][0] * x[4][0] - x[5][0] - sinh(2. * x[7][0])* x[8][0] + 2. * x[9][0] + 2.0004339741653854440,

		sin(x[0][0] * x[1][0]) + x[2][0] * x[8][0] * x[6][0] - exp(-x[9][0] + x[5][0]) + 3. * x[4][0] * x[4][0] - x[5][0] * (x[7][0] + 1) + 10.886272036407019994,

		x[0][0] - x[1][0] + x[2][0] - x[3][0] + x[4][0] - x[5][0] + x[6][0] - x[7][0] + x[8][0] - x[9][0] - 3.1361904761904761904,

		2. * cos(-x[8][0] + x[3][0]) + x[4][0] / (x[2][0] + x[0][0]) - sin(x[1][0] * x[1][0]) + cos(x[6][0] * x[9][0])*cos(x[6][0] * x[9][0]) - x[7][0] - 0.1707472705022304757,

		sin(x[4][0]) + 2. * x[7][0] * (x[2][0] + x[0][0]) - exp(-x[6][0] * (-x[9][0] + x[5][0])) + 2. * cos(x[1][0]) - 1. / (x[3][0] - x[8][0]) - 0.3685896273101277862,

		exp(x[0][0] - x[3][0] - x[8][0]) + x[4][0] * x[4][0] / x[7][0] + 0.5*cos(3 * x[9][0] * x[1][0]) - x[5][0] * x[2][0] + 2.0491086016771875115,

		x[1][0] * x[1][0] * x[1][0] * x[6][0] - sin(x[9][0] / x[4][0] + x[7][0]) + (x[0][0] - x[5][0]) * cos(x[3][0]) + x[2][0] - 0.7380430076202798014,

		x[4][0] * (x[0][0] - 2. * x[5][0])*(x[0][0] - 2. * x[5][0]) - 2. * sin(-x[8][0] + x[2][0]) + 1.5*x[3][0] - exp(x[1][0] * x[6][0] + x[9][0]) + 3.5668321989693809040,

		7. / x[5][0] + exp(x[4][0] + x[3][0]) - 2 * x[1][0] * x[7][0] * x[9][0] * x[6][0] + 3. * x[8][0] - 3. * x[0][0] - 8.4394734508383257499,

		x[9][0] * x[0][0] + x[8][0] * x[1][0] - x[7][0] * x[2][0] + sin(x[3][0] + x[4][0] + x[5][0]) *x[6][0] - 0.7823809523809523809,

	};
	return result=SNAE_fm;
};

auto jacobi_matrix = [](const Matrix &x)
{
	Matrix result(10, 10);
	
	double Jac_form[] = {
		-sin(x[0][0] * x[1][0]) * x[1][0],
		-sin(x[0][0] * x[1][0]) * x[0][0],
		3. * exp(-(3 * x[2][0])),
		x[4][0] * x[4][0],
		2 * x[3][0] * x[4][0],
		-1,
		0,
		-2. * cosh((2 * x[7][0])) * x[8][0],
		-sinh((2 * x[7][0])),
		2,
		cos(x[0][0] * x[1][0]) * x[1][0],
		cos(x[0][0] * x[1][0]) * x[0][0],
		x[8][0] * x[6][0],
		0,
		6 * x[4][0],
		-exp(-x[9][0] + x[5][0]) - x[7][0] - 0.1e1,
		x[2][0] * x[8][0],
		-x[5][0],
		x[2][0] * x[6][0],
		exp(-x[9][0] + x[5][0]),
		1,
		-1,
		1,
		-1,
		1,
		-1,
		1,
		-1,
		1,
		-1,
		-x[4][0] * pow(x[2][0] + x[0][0], -2.),
		-2. * cos(x[1][0] * x[1][0]) * x[1][0],
		-x[4][0] * pow(x[2][0] + x[0][0], -2.),
		-2. * sin(-x[8][0] + x[3][0]),
		1. / (x[2][0] + x[0][0]),
		0,
		-2. * cos(x[6][0] * x[9][0]) * sin(x[6][0] * x[9][0]) * x[9][0],
		-1,
		2. * sin(-x[8][0] + x[3][0]),
		-2. * cos(x[6][0] * x[9][0]) * sin(x[6][0] * x[9][0]) * x[6][0],
		2 * x[7][0],
		-2. * sin(x[1][0]),
		2 * x[7][0],
		pow(-x[8][0] + x[3][0], -2.),
		cos(x[4][0]),
		x[6][0] * exp(-x[6][0] * (-x[9][0] + x[5][0])),
		-(x[9][0] - x[5][0]) * exp(-x[6][0] * (-x[9][0] + x[5][0])),
		(2 * x[2][0]) + 2. * x[0][0],
		-pow(-x[8][0] + x[3][0], -2.),
		-x[6][0] * exp(-x[6][0] * (-x[9][0] + x[5][0])),
		exp(x[0][0] - x[3][0] - x[8][0]),
		-3. / 2. * sin(3. * x[9][0] * x[1][0]) * x[9][0],
		-x[5][0],
		-exp(x[0][0] - x[3][0] - x[8][0]),
		2 * x[4][0] / x[7][0],
		-x[2][0],
		0,
		-x[4][0] * x[4][0] * pow(x[7][0], (-2)),
		-exp(x[0][0] - x[3][0] - x[8][0]),
		-3. / 2. * sin(3. * x[9][0] * x[1][0]) * x[1][0],
		cos(x[3][0]),
		3. * x[1][0] * x[1][0] * x[6][0],
		1,
		-(x[0][0] - x[5][0]) * sin(x[3][0]),
		cos(x[9][0] / x[4][0] + x[7][0]) * x[9][0] * pow(x[4][0], (-2)),
		-cos(x[3][0]),
		pow(x[1][0], 3.),
		-cos(x[9][0] / x[4][0] + x[7][0]),
		0,
		-cos(x[9][0] / x[4][0] + x[7][0]) / x[4][0],
		2. * x[4][0] * (x[0][0] - 2. * x[5][0]),
		-x[6][0] * exp(x[1][0] * x[6][0] + x[9][0]),
		-2. * cos(-x[8][0] + x[2][0]),
		0.15e1,
		pow(x[0][0] - 2. * x[5][0], 2.),
		-4. * x[4][0] * (x[0][0] - 2. * x[5][0]),
		-x[1][0] * exp(x[1][0] * x[6][0] + x[9][0]),
		0,
		2. * cos(-x[8][0] + x[2][0]),
		-exp(x[1][0] * x[6][0] + x[9][0]),
		-3,
		-2. * x[7][0] * x[9][0] * x[6][0],
		0,
		exp((x[4][0] + x[3][0])),
		exp((x[4][0] + x[3][0])),
		-0.7e1 * pow(x[5][0], -2.),
		-2. * x[1][0] * x[7][0] * x[9][0],
		-2. * x[1][0] * x[9][0] * x[6][0],
		3,
		-2. * x[1][0] * x[7][0] * x[6][0],
		x[9][0],
		x[8][0],
		-x[7][0],
		cos(x[3][0] + x[4][0] + x[5][0]) * x[6][0],
		cos(x[3][0] + x[4][0] + x[5][0]) * x[6][0],
		cos(x[3][0] + x[4][0] + x[5][0]) * x[6][0],
		sin(x[3][0] + x[4][0] + x[5][0]),
		-x[2][0],
		x[1][0],
		x[0][0]
		

	};
	
	return result=Jac_form;
};

double Pol(double x)
{
	return tan(0.5*x + 0.2) - x*x;
}
double Df(double(*Func)(double), double x)
{
	return (Func(x + E) - Func(x)) / E;
}
double Newton(double(*Func)(double), double x0)
{

	double x1=x0-Func(x0)/ Df(Func, x0);

	int Count = 0;
	while (abs(x1 - x0) > Eps && Count < 200)
	{
		
		++Count;		
		x0 = x1;
		x1 = x1 - Func(x1) / Df(Func, x1);
	}
	std::cout << "iter " << Count<<std::endl;
	return x1;
}
bool checkMatrix(const Matrix &A, const Matrix &B)
{
	for (int i = 0; i < A.dimM; ++i)
	{
		if ((A[i][0] - B[i][0]) > Eps) return true;
	}
	return false;
};
Matrix NewtonForSystem(Matrix X0, int ModifIterationsNum=-1, bool Hybrid=false)
{
	bool Modif = false;
	Matrix X1(10, 1);
	Matrix JacobiX0(10, 10);
	JacobiX0 = inverse(jacobi_matrix(X0));

	X1=X0- JacobiX0*SNAE_matrix(X0);
		
	

	int operations = 3;
	int Count = 0;
	while (checkMatrix(X1,X0)&& Count < 500)
	{
		
		if (Hybrid&&ModifIterationsNum > 0)
			//------------------------Hybrid----------------------------------
		{
			Modif = true;
			if (Count%ModifIterationsNum == 0 && Count > 0)
			{
				JacobiX0 = inverse(jacobi_matrix(X1));
				operations += 1;
			}
		
		}
		//---------------------------------------------------------------
		else
			//------------------------Modif----------------------------------
			if (Count == ModifIterationsNum)
				Modif = true;
		//---------------------------------------------------------------
	
		++Count;
		X0 = X1;
		++operations;
			
		/*	std::cout << "Jac\n";
			jacobi_matrix(X1).show();
			std::cout << "inverse\n";
			(inverse(jacobi_matrix(X1))).show();
			(inverse(jacobi_matrix(X1))*jacobi_matrix(X1)).show();
			std::cout << "SNAE\n";
			SNAE_matrix(X1).show();*/
		if (!Modif)
		{
			X1 = X1 - inverse(jacobi_matrix(X1))*SNAE_matrix(X1);
			operations += 4;
		}
	
		else
		{
			X1 = X1 - JacobiX0*SNAE_matrix(X1);
			operations += 3;
		}
		
		//	X1.show();
			//SNAE_matrix(X1).show();
	}
		std::cout << "Count: " << Count<<std::endl;
		std::cout << "Operations: " << operations<< std::endl;
		return X1;
}
void getInfForSolve(Matrix &A,int iter, bool Hybrid=false)
{
	Matrix ans(10, 1);

	clock_t start = clock();
	ans = NewtonForSystem(A,iter,Hybrid);
	clock_t end = clock();
	double seconds = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "time: " << seconds << std::endl;
	ans.show();
	//(SNAE_matrix(ans)).show();
	if (abs(SNAE_matrix(ans)[0][0])<Eps) std::cout << "ans correct \n";
	else  std::cout << "ans uncorrect! \n";
	
}
int main()
{
	std::cout.precision(17);
	int x0 = 1;
	std::cout << Newton(Pol, x0) << std::endl;
	std::cout << Pol(Newton(Pol, x0)) << std::endl;
	std::cout << "---------------------------------------------\n";
	

	for (int i = -1; i <8; i++)
	{
		std::cout << i << "\n";
		getInfForSolve(Matrix_x0(),i);
		std::cout << std::endl << "Hybrid: \n";
		getInfForSolve(Matrix_x0(), i, true);
		std::cout << std::endl;
	}

	
	/*std::cout << "\n\n\n" << "for X02: \n";
	std::cout << "---------------------------------------------\n";

	for (int i = -1; i <5; i++)
	{
		std::cout << i << "\n";
		getInfForSolve(Matrix_x02(), i);
		std::cout << std::endl<<"Hybrid: \n";
		getInfForSolve(Matrix_x02(), i, true);
		std::cout << std::endl;
	}
	
	return 0;*/
	

}