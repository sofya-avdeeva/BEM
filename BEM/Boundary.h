#pragma once
#define _USE_MATH_DEFINES

#include <vector>
#include <string>
#include <fstream>
#include <functional>

#include "BEM.h"
#include "Math.h"
#include "Point.h"

using namespace std;

struct Edge
{
	int v1, v2;
	int valueNo;
};

void InputBound(vector<Edge>& bound, string boundsFile);
void Boundary1(vector<vector<double>>& S, vector<double>& b, vector<Edge>& bound, vector<Point>& points, bool hasOffset);
void Boundary2(vector<double>& b, vector<Edge>& bound, vector<Point>& points, bool hasOffset);