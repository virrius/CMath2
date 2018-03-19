#pragma once
#include<iostream>
#include"Matrix.h"
#include"func.h"
#include<string>
#include"tests.h"



int main()
{
	
	Matrix  A(3,3),PR(A.dimM, A.dimM), PC(A.dimM, A.dimM), L(A.dimM, A.dimM), U(A.dimM, A.dimM), B(A.dimM, 1), X(A.dimM, 1);
	std::cout << testPLU(4, 4, 1)<< testPLU(5, 5, 0)<<testPLU(3,3,0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testSolve(3,1) << testSolve(4,0)<< testSolve(8, 0);
	std::cout << "----------------------------------------------------\n";
	std::cout << testInv(3,3,1)<<testInv(5,5,0)<<testInv(8,8,0);
	A.Ex4();

	std::cout<<A.Norm()<<" "<<inverse(A).Norm()<<" "<<A.CondNum();

	return 0;
}