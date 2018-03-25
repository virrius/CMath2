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
		B.show();


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
std::string testSolveQR(int dimM, int dimN, bool debug)
{
	Matrix A(dimM, dimN), B(A.dimM, 1), X(A.dimM, 1);
		
	A.random();
	B.random();

		SolveQR(A, B, X);
	if (debug)
	{
		std::cout << " matrix A\n";
		A.show();
	
		std::cout << " matrix x\n";
		X.show();
		std::cout << " matrix AX \n";
		(A*X).show();
		std::cout << " matrix B\n";
		B.show();
		
	}
	if (A*X==B)
		return "SolveQR ok\n";
	else
		return "SolveQR wrong\n";
};
std::string testQR(int dimM, int dimN, bool debug)
{
	Matrix A(dimM, dimN), Q(dimM, dimN), R(dimM, dimN);

	A.random();


	QR(A, Q, R);
	if (debug)
	{
		std::cout << " matrix A\n";
		A.show();


		std::cout << " matrix Q \n";
		Q.show();
		std::cout << " matrix R \n";
		R.show();
		std::cout << " matrix QR \n";
		(Q*R).show();
		std::cout << " matrix A \n";
		A.show();

	}
	if (A == Q*R)
		return "QR ok\n";
	else
		return "QR wrong\n";
};
std::string testJacobi(int dimM, int dimN, bool debug)
{
	Matrix A(dimM, dimN), B(A.dimM, 1), X(A.dimM, 1);

	A.randomDP();
	B.random();

	SolveJacobi(A, B, X);
	if (debug)
	{
		std::cout << " matrix A\n";
		A.show();

		std::cout << " matrix x\n";
		X.show();
		std::cout << " matrix AX \n";
		(A*X).show();
		std::cout << " matrix B\n";
		B.show();

	}
	if (A*X == B)
		return "SolveJacobi ok\n";
	else
		return "SolveJacobi wrong\n";
};
std::string testZ(int dimM, int dimN, bool debug)
{
	Matrix A(dimM, dimN), B(A.dimM, 1), X(A.dimM, 1);

	A.randomDP();
	B.random();

	SolveZ(A, B, X);
	if (debug)
	{
		std::cout << " matrix A\n";
		A.show();

		std::cout << " matrix x\n";
		X.show();
		std::cout << " matrix AX \n";
		(A*X).show();
		std::cout << " matrix B\n";
		B.show();

	}
	if (A*X == B)
		return "SolveZ ok\n";
	else
		return "SolveZ wrong\n";
};