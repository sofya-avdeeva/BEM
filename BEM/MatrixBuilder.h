#pragma once
#define _USE_MATH_DEFINES

#include <math.h>
#include <vector>

#include "MeshBuilder.h"
#include "BEM.h"
#include "Math.h"
#include "mkl.h"

using namespace std;

void Build(
	vector<vector<double>>& V, vector<vector<double>>& K, vector<vector<double>>& Kt,
	vector<vector<double>>& D, vector<Point>& points, vector<BEMElement>& elements);
void ShurComplement(
	vector<vector<double>>& V, vector<vector<double>>& K,
	vector<vector<double>>& D, vector<vector<double>>& S,
	int phiCount, int psiCount);

double* MatrixToVector(int n, int m, vector<vector<double>> vec);
vector<vector<double>> VectorToMatrix(int n, int m, double* matrix);
