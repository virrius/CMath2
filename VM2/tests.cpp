#include"tests.h"


std::string testInv(int dimM, int dimN, bool debug)
{

	Matrix A(dimM, dimN), PR(A.dimM, A.dimM), PC(A.dimM, A.dimM), L(A.dimM, A.dimM), U(A.dimM, A.dimM), B(A.dimM, 1), X(A.dimM, 1);
	A.random();
	Matrix inv = inverse(A);
	Matrix idn(A.dimM, A.dimN);
	idn.Identity();

	if (debug)
	{
		std::cout << " matrix A\n";
		A.show();
		std::cout << " matrix inv\n";
		inv.show();
		std::cout << " matrix A*inv\n";

		PLU(A, L, U, PC, PR);
		(A*inv).show();

	}
	if (A*inv == idn)
		return "Inv ok\n";
	else
		return "Inv wrong\n";
}
std::string testSolve(int dim, bool debug)
{
	Matrix A(dim, dim), PR(dim, dim), PC(dim, dim), L(dim, dim), U(dim, dim), B(dim, 1), X(dim, 1);
	A.random();
	B.random();
	PLU(A, L, U, PC, PR);
	Solve(L, U, B, X);
	if (debug)
	{
		std::cout << " Matrix A";
		A.show();

		std::cout << " Matrix X";
		X.show();
		std::cout << " Matrix A*X";
		(PR*A*PC*X).show();
		std::cout << " Matrix B";
		(PR.Transpose()*B).show();


	}
	if (PR*A*PC*X == B)
		return "Solve ok\n";
	else
		return "Solve wrong\n";
}
std::string testPLU(int dimM, int dimN, bool debug)
{
	Matrix A(dimM, dimN), PR(dimM, dimN), PC(dimM, dimN), L(dimM, dimN), U(dimM, dimN);

	A.random();


	PLU(A, L, U, PC, PR);
	if (debug)
	{
		std::cout << " matrix A\n";
		A.show();
		std::cout << " matrixes P\n";
		PR.show();
		PC.show();

		std::cout << " matrix L \n";
		L.show();
		std::cout << " matrix U \n";
		U.show();
		std::cout << " matrix PLU \n";
		(PR.Transpose()*L*U*PC.Transpose()).show();
		std::cout << " matrix A \n";
		A.show();

	}
	if (L*U == PR*A*PC)
		return "PLU ok\n";
	else
		return "PLU wrong\n";
};