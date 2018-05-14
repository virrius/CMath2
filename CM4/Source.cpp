#include<iostream>
#include"Matrix.h"
#include"func.h"
double nu0()
{
	return 5.8955591631161391715427427670185;

}
double nu1()
{
	return 12.212229695026288283909967160253;

}
double nu2()
{
	return 25.56664299495079122270233602361;
}
Matrix Xsj(double &x1, double & x2, double & x3)
{
	double m[] = { 1,1,1,x1,x2,x3,x1*x1,x2*x2,x3*x3};
	Matrix B(3, 3);
	return B = m;
}
Matrix Nu()
{
	double n[] = { nu0(),nu1(),nu2() };
	Matrix Nu(3, 1);
	return Nu = n;
}
double func(double &x)
{
	return  4*cos(0.5*x)* exp(-5*x / 4) + 2*sin(4.5*x)* exp(x / 8) + 2;
}
double weightFunc(double &x)
{
	double a = 1.3;
	double b = 2.2;
	return pow(x - a, 0)*pow(b - x, 5. / 6);
}
void intQuadForm()
{
	double x1 = 1.3;
	double x2 = (1.3 + 2.2) / 2;
	double x3 = 2.2;
	auto A = Xsj(x1, x2, x3);
	auto B = Nu();
	Matrix X(3, 1);
	SolveQR(A, B, X);
	double integral = X[0][0]*func(x1)+X[1][0]*func(x2)+
}


int main()
{
	intQuadForm();


	std::cin.get();
	return 0;
}