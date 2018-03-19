#pragma once

#include<iostream>
#include"Matrix.h"





void FindLeadElement(Matrix &A, Matrix &PC, Matrix &PR);
void Solve(Matrix &L, Matrix &U, Matrix &B, Matrix &X);
void PLU(Matrix &A, Matrix &L, Matrix &U, Matrix &PC, Matrix &PR);
Matrix inverse(Matrix &A);