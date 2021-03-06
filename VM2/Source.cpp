#pragma once
#include<iostream>
#include"Matrix.h"
#include"func.h"
#include<string>
#include"tests.h"
#include<time.h>


int main()
{
	//srand(time(0));
	Matrix  A(3, 3), PR(A.dimM, A.dimM), PC(A.dimM, A.dimM), L(A.dimM, A.dimM), U(A.dimM, A.dimM), B(A.dimM, 1), X(A.dimM, 1), R(A.dimM, A.dimN), Q(A.dimM, A.dimN);
	std::cout << testPLU(4, 4, 1)<< testPLU(5, 5, 0)<<testPLU(3,3,0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testSolvePLU(3,1) << testSolvePLU(4,0)<< testSolvePLU(8, 0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testInv(3,3,1)<<testInv(5,5,0)<<testInv(8,8,0);
	std::cout << "----------------------------------------------------\n";
	A.Ex4();
	std::cout<<A.Norm()<<" "<<inverse(A).Norm()<<" "<<A.CondNum()<<std::endl;
	std::cout << "----------------------------------------------------\n";
	std::cout << testQR(3, 3, 1) << testQR(5, 5, 0) << testQR(8, 8, 0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testSolveQR(4, 4, 1) << testSolveQR(6, 6, 0) << testSolveQR(8,8,0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testJacobi(8, 8, 1)<<testJacobi(4,4,0)<<testJacobi(6,6,0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testZ(4, 4, 1)<<testZ(6,6,0)<<testZ(8,8,0);
	std::cout << "----------------------------------------------------\n";
	JacobiVsZ();
	A.EXJac();
	B.random(); 
	int i;
	SolveJacobi(A, B, X, i);
	X.show();
	return 0;
}