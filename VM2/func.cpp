#include"func.h";
#include"Matrix.h"

void FindLeadElement(Matrix &A, Matrix &PC, Matrix &PR)
{
for (int i = 0; i < A.dimM; i++)
{
	double max = 0;
	int  swapR = -1, swapC = -1;
	for (int row = i; row < A.dimM; row++)
	{
		for (int col = i; col < A.dimN; col++)
			if (std::abs(A[row][col]) > std::abs(max))
			{
				max = A[row][col];
				swapR = row;
				swapC = col;
			}
		//std::cout <<"swapR "<< swapR << " " <<swapC << " " << i << std::endl;

	}
	if (max == 0)
	{
		std::cout << "Null";
		return;
	}
	PR.SwapRows(swapR, i);
	A.SwapRows(swapR, i);

	PC.SwapCol(swapC, i);
	A.SwapCol(swapC, i);
}
/*std::cout << "after swap ";
A.show();
P.show();*/

};
void Solve(Matrix &L, Matrix &U, Matrix &B, Matrix &X)
{
	int n = B.dimM - 1;
	Matrix Y(B.dimM, B.dimN);

	for (int i = 0; i <= n; i++)
	{
		double sum = 0;
		for (int k = 0; k < i; k++) {
			sum += L[i][k] * Y[k][0];

		}
		Y[i][0] = B[i][0] - sum;

	}

	//Y.show();
	for (int i = 0; i <= n; i++)
	{
		double sum = 0;
		for (int k = 0; k < i; k++) {
			sum += U[n - i][n - k] * X[n - k][0];
		}
		X[n - i][0] = (Y[n - i][0] - sum) / U[n - i][n - i];

	}



}
void PLU(Matrix &A, Matrix &L, Matrix &U, Matrix &PC, Matrix &PR)
{
	PC.Identity();
	PR.Identity();
	U = A;




	FindLeadElement(U, PC, PR);

	for (int k = 0; k < A.dimN; k++)
	{

		for (int i = k; i < A.dimN; i++)
		{

			for (int j = i; j < A.dimM; j++)
			{

				L[j][i] = U[j][i] / U[i][i];
			}
		}
		//std::cout << "matrix L " << k;
		//L.show();

		for (int i = k + 1; i < A.dimM; i++)
			for (int j = k; j < A.dimN; j++)
				U[i][j] = U[i][j] - L[i][k] * U[k][j];
		//std::cout << "matrix U " << k;
		//U.show();
	}
}
Matrix inverse(Matrix &A)
{
	Matrix  PR(A.dimM, A.dimM), PC(A.dimM, A.dimM), L(A.dimM, A.dimM), U(A.dimM, A.dimM), B(A.dimM, 1), X(A.dimM, 1);

	Matrix Res(A.dimM, A.dimN), Tmp(A.dimM, 1);
	PLU(A, L, U, PC, PR);
	//PR.show();
	//	PC.show();

	for (int i = 0; i < A.dimM; i++)
	{
		if (i != 0)B[i - 1][0] = 0;
		B[i][0] = 1;
		Solve(L, U, B, Tmp);
		//std::cout << "TMP\n";
		//Tmp.show();
		for (int j = 0; j < A.dimM; j++)
		{
			Res[j][i] = Tmp[j][0];
		}
	}


	return PC*Res*PR;
};