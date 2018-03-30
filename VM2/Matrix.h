#pragma once


const double E = 0.00001;
				
class Matrix {
	double **mas;
	int sign = 1;
public:
	int dimM,dimN;

	//Matrix( Matrix &A)
	//{
	//	int n = std::fmax(A.dimM, A.dimN);
	//	dimM = n;
	//	dimN = n;
	//	mas = new double*[n];
	//	for (int i = 0; i < n; i++)
	//		mas[i] = new double[n];
	//	for (int i = 0; i < n; i++)
	//		for (int j = 0; j < n; j++)
	//			mas[i][j] = 0;
	//}
	Matrix(int M, int N);
	//Matrix::~Matrix();
	void show()const;
	double* operator[](int i)const;
	bool Matrix::operator==(const Matrix &A);
	
	const Matrix operator=(const Matrix &A);
	
	Matrix operator*(const Matrix &A);
	Matrix operator*(double K);
	Matrix operator+(const Matrix &A);
	void randomDP();
	void random() {
		for (int i = 0; i < dimM; i++)
			for (int j = 0; j < dimN; j++)
				mas[i][j] = rand() % 100 - 50;
	}
	void Identity() {
		for (int i = 0; i < dimM; i++)
			for (int j = 0; j < dimN; j++)
				if (i == j)
					mas[i][j] = 1;
				else
				{
					mas[i][j] = 0;
				}
	}
	void SwapRows(int f, int s);	
	void SwapCol(int f, int s);
	
	Matrix Transpose();
	double Determinant();
	int Rank(Matrix &U);	
	double CondNum();
	double Norm();
	//Matrix examples;
	void Ex1()
	{
		mas[0][0] = 10;
		mas[0][1] = -7;
		mas[0][2] = 0;
		mas[1][0] = -3;
		mas[1][1] = 6;
		mas[1][2] = 2;
		mas[2][0] = 5;
		mas[2][1] = -1;
		mas[2][2] = 5;
	}
	void Ex2()
	{
		mas[0][0] = 8;
		mas[0][1] = 4;
		mas[0][2] = 2;
		mas[0][3] = 1;
		mas[1][0] = 16;
		mas[1][1] = 9;
		mas[1][2] = 3;
		mas[1][3] = 1;
		mas[2][0] = 32;
		mas[2][1] = 93;
		mas[2][2] = 42;
		mas[2][3] = 8;
		mas[3][0] = 48;
		mas[3][1] = 5;
		mas[3][2] = 12;
		mas[3][3] = 2;



	}
	void Ex3()
	{
		mas[0][0] = 2;
		mas[0][1] = 0;
		mas[0][2] = 2;
		mas[0][3] = 0.6;
		mas[1][0] = 3;
		mas[1][1] = 3;
		mas[1][2] = 4;
		mas[1][3] = -2;
		mas[2][0] = 5;
		mas[2][1] = 5;
		mas[2][2] = 4;
		mas[2][3] = 2;
		mas[3][0] = -1;
		mas[3][1] = -2;
		mas[3][2] = 3.4;
		mas[3][3] = -1;


	}
	void Ex4()
	{
		mas[0][0] = 2;
		mas[0][1] = 7;
		mas[0][2] = 6;
		mas[1][0] = 9;
		mas[1][1] = 5;
		mas[1][2] = 1;
		mas[2][0] = 4;
		mas[2][1] = 3;
		mas[2][2] = 8;

	}	
	void ExSLAEA()
	{

		mas[0][0] = 8;
		mas[0][1] = 1;
		mas[0][2] = 1;
		mas[1][0] = 1;
		mas[1][1] = 10;
		mas[1][2] = 1;
		mas[2][0] = 1;
		mas[2][1] = 1;
		mas[2][2] = 12;
	
	}
	void ExSLAEB()
	{

		mas[0][0] = 10;
		mas[1][0] = 12;
		mas[2][0] = 14;
	}
	void EXJac()
	{
		mas[0][0] = 2;
		mas[0][1] = -1;
		mas[0][2] = 2;
		mas[1][0] = -1;
		mas[1][1] = 1;
		mas[1][2] = -3;
		mas[2][0] = 2;
		mas[2][1] = -3;
		mas[2][2] = 11;
	}
	
	
};




