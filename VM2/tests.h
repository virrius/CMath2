#pragma once
#include<string>
#include"Matrix.h"
#include"func.h"

std::string testInv(int dimM, int dimN, bool debug);
std::string testSolve(int dim, bool debug);
std::string testPLU(int dimM, int dimN, bool debug);
std::string testQR(int dimM, int dimN, bool debug);
std::string testSolveQR(int dimM, int dimN, bool debug);
std::string testJacobi(int dimM, int dimN, bool debug);
std::string testZ(int dimM, int dimN, bool debug);