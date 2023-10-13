#pragma once

#include <vector>
#include <functional>

#include "Point.h"

using namespace std;

const vector<function<double(double)>> psi = {
	[](double r) { return 1.0 - r; },
	[](double r) { return r; }
};

const vector<function<double()>> phi = {
	[]() { return 1.0; }
};

const vector<function<double()>> gradPsi = {
	[]() { return -1; },
	[]() { return 1; }
};

const vector<function<double(double, double)>> uValueQ = {
	[](double r, double z) { return z; },
	[](double r, double z) { return z; }
};

const vector<function<double(double, double)>> thetaValueQ = {
	[](double r, double z) { return 0.0; },
	[](double r, double z) { return 0.0; }
};

inline double uQ(double r, double z)
{
	return 0;
}