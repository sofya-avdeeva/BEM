#include "Boundary.h"

void InputBound(vector<Edge>& bound, string boundsFile)
{
	int boundCount, numCond;
	int edgeCount;
	double beg, end;
	int n;
	ifstream in(boundsFile);
	in >> boundCount;

	for (int i = 0; i < boundCount; i++)
	{
		in >> beg >> end >> n >> numCond;

		for (int j = 0; j < n; j++)
		{
			Edge b;
			b.v1 = beg + j;
			b.v2 = beg + 1 + j;
			b.valueNo = numCond;

			bound.push_back(b);
		}
	}

	if (end == 0)
		bound[bound.size() - 1].v2 = 0;

	in.close();
}

void Boundary1(vector<vector<double>>& S, vector<double>& b, vector<Edge>& bound, vector<Point>& points, bool hasOffset)
{
	int offset = hasOffset ? points.size() : 0;

	for (auto& edge : bound)
	{
		function<double(double, double)> Ug = uValueQ[edge.valueNo];
		vector<double> f(2);

		double r0 = points[edge.v1].r;
		double r1 = points[edge.v2].r;
		double z0 = points[edge.v1].z;
		double z1 = points[edge.v2].z;

		double h = Distance(points[edge.v2], points[edge.v1]);

		double a1 = r0;
		double b1 = r1;

		if (a1 == b1)
		{
			a1 = z0;
			b1 = z1;
		}

		bool needReverse = a1 > b1;

		if (needReverse)
			swap(a1, b1);

		vector<vector<double>> M(2, vector<double>(2));

		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				M[i][j] += h * Gauss18([&](double ksi) { return psi[i](ksi) * psi[j](ksi) * (r0 + ksi * (r1 - r0)); });

		for (int i = 0; i < 2; i++)
			f[i] = Gauss18([&](double ksi) { return Ug(r0 + ksi * (r1 - r0), z0 + ksi * (z1 - z0)) * psi[i](ksi); });

		vector<double> q(2);
		Gauss(M, q, f);

		S[edge.v1 + offset][edge.v1 + offset] = 1.0E+50;
		b[edge.v1 + offset] = 1.0E+50 * Ug(r0, z0);
		S[edge.v2 + offset][edge.v2 + offset] = 1.0E+50;
		b[edge.v2 + offset] = 1.0E+50 * Ug(r1, z1);
	}
}

void Boundary2(vector<double>& b, vector<Edge>& bound, vector<Point>& points, bool hasOffset)
{
	int offset = hasOffset ? points.size() : 0;

	for (auto& edge : bound)
	{
		vector<int> v = { edge.v1, edge.v2 };

		function<double(double, double)> theta = thetaValueQ[edge.valueNo];
		double h = Distance(points[v[0]], points[v[1]]);

		double r0 = points[v[0]].r;
		double r1 = points[v[1]].r;
		double z0 = points[v[0]].z;
		double z1 = points[v[1]].z;

		vector<double> a(2);

		for (int i = 0; i < 2; i++)
			b[v[i] + offset] += h * Gauss18([&](double ksi) { return theta(r0 + ksi * (r1 - r0), z0 + ksi * (z1 - z0)) * psi[i](ksi) * (r0 + ksi * (r1 - r0)); });

		double tem = 0.0;
	}
}