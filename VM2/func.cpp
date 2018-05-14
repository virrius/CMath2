#pragma once
#include"func.h"
#include"Matrix.h"
void FindLeadElement(Matrix &A, Matrix &PC, Matrix &PR, int Minor)
{
for (int i = Minor; i < A.dimM; i++)
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
void SolvePLU(Matrix &L, Matrix &U, Matrix &B, Matrix &X)
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




	

	for (int k = 0; k < A.dimN; k++)
	{
		FindLeadElement(U, PC, PR, k);
		//std::cout << "matrix L " << k;
		//L.show();
		//U.show();
		
		for (int i = k+1 ; i < A.dimM; i++)
		{
			U[i][k] = U[i][k] / U[k][k];
			for (int j = k+1; j < A.dimN; j++)
				U[i][j] = U[i][j] - U[i][k] * U[k][j];
		}
		
		//std::cout << "matrix U " << k;
		//U.show();
	}
	for (int i = 0; i < A.dimM; ++i)
	{
		L[i][i] = 1;
		for (int j = 0; j < A.dimM; ++j)
			if (i > j)
			{
				L[i][j] = U[i][j];
				U[i][j] = 0;
			}
			
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
		SolvePLU(L, U, B, Tmp);
		//std::cout << "TMP\n";
		//Tmp.show();
		for (int j = 0; j < A.dimM; j++)
		{
			Res[j][i] = Tmp[j][0];
		}
	}


	return PC*Res*PR;
};
void QR(Matrix &A, Matrix &Q, Matrix &R)
{
	Matrix T(A.dimM, A.dimM);
	Q.Identity();
	R = A;
	for (int i = 0; i < A.dimM - 1; i++)
	{

		for (int j = i + 1; j < A.dimM; j++)
		{
			if (R[j][i] != 0)
			{
				T.Identity();
				T[i][i] = R[i][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				T[i][j] = R[j][i] / sqrt(R[i][i] * R[i][i] + R[j][i] * R[j][i]);
				T[j][j] = T[i][i];
				T[j][i] = -T[i][j];

			}
			R = T*R;
			Q = T*Q;
		//	T.show();
			//R.show();
		}


	}
	Q = inverse(Q);

}
void SolveQR(Matrix &A, Matrix &B, Matrix &X)
{
	Matrix Q(A.dimM, A.dimN), R(A.dimM, A.dimN),Y(X.dimM,X.dimN);
	QR(A, Q, R);
	Q = inverse(Q);
	Y = Q*B;
	//Y.show();
// R.show();
	X = inverse(R)*Y;
	//X.show();
}
void SolveJacobi(Matrix &A, Matrix &B, Matrix &X,int &Debug_Count)
{
	const double Eps =  0.00000001;
	Matrix D(A.dimM, A.dimN), invD(D.dimM, D.dimN), G(D.dimM, 1), C(A.dimM, A.dimN),newX(X.dimM,X.dimN);

	for (int i = 0; i < A.dimM; i++)
	{
		D[i][i] = A[i][i];
		invD[i][i] = 1 / D[i][i];
		
	}
	
	C = invD*(D + A*(-1));
	G = invD*B;
	X = G;
	Debug_Count = 0;
	double norm;
	do {
		norm = 0;

		newX = C*X + G;
	
		for (int i = 0; i < X.dimM; i++)
		{
			
			if (abs(X[i][0]- newX[i][0]) > norm)
				norm = abs(X[i][0] - newX[i][0]);
		}
		X = newX;
		//X.show();

		Debug_Count++;
		if (Debug_Count == 1000)
		{
			std::cout<< "Jacobi error \n";
			return;
		}
	} while (norm>Eps);	
	
}
void SolveZ(Matrix &A, Matrix &B, Matrix &X, int &Debug_Count)
{
	const double Eps = 0.00000001;
	Matrix D(A.dimM, A.dimN), U(D.dimM, D.dimN), L(D.dimM, D.dimN), invLD(A.dimM, A.dimN), newX(X.dimM, X.dimN);
	for (int i = 0; i < A.dimM; i++)
		for (int j = 0; j < A.dimM; j++)
		{
			if (i == j) D[i][j] = A[i][j];
			if (i < j) U[i][j] = A[i][j];
			if (i > j)L[i][j] = A[i][j];
		}
	invLD = inverse(L + D);

	double norm;
	Debug_Count = 0;
	do {
		norm = 0;
		newX = invLD*(U*(-1)*X + B);
		//newX.show();
		for (int i = 0; i < A.dimM; i++)
		{
			norm += (newX[i][0] - X[i][0])*(newX[i][0] - X[i][0]);
		}
		norm = sqrt(norm);

		X = newX;
		Debug_Count++;
	} while (norm > Eps);
}
