#pragma once

#include<iostream>
#include"Matrix.h"




void QR(Matrix &A, Matrix &Q, Matrix &R);
void FindLeadElement(Matrix &A, Matrix &PC, Matrix &PR, int Minor);
void SolvePLU(Matrix &L, Matrix &U, Matrix &B, Matrix &X);
void PLU(Matrix &A, Matrix &L, Matrix &U, Matrix &PC, Matrix &PR);
Matrix inverse(Matrix &A);
void SolveQR(Matrix &A, Matrix &B, Matrix &X);
void SolveJacobi(Matrix &A, Matrix &B, Matrix &X, int &Debug_Count);
void SolveZ(Matrix &A, Matrix &B, Matrix &X, int &Debug_Count);
