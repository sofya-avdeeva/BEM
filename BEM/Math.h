#pragma once

#include <vector>

#include "Point.h"
#include <iostream>

using namespace std;

inline Point Cross(Point a)
{
	return Point(a.z, -a.r);
}

inline Point Normal(Point a, Point b)
{
	return Cross(b - a).Normalize();
}

inline double Distance(Point a, Point b)
{
	return sqrt((a.r - b.r) * (a.r - b.r) + (a.z - b.z) * (a.z - b.z));
}

struct QuadratureNode
{
   QuadratureNode(Point node, double weight) : Node(node), Weight(weight) {}
   Point Node;
   double Weight;
};

inline vector<QuadratureNode> SegmentGaussOrder23()
{
	const int n = 12;
	double x[n] = { -0.981560634246719, -0.904117256370475, -0.769902674194305, -0.587317954286617, 
		-0.367831498998180, -0.125233408511469, 0.125233408511469, 0.367831498998180,
		0.587317954286617, 0.769902674194305, 0.904117256370475, 0.981560634246719 };
	double w[n] = { 0.047175336386512, 0.106939325995318, 0.160078328543346, 0.203167426723066,
		0.233492536538355, 0.249147045813403, 0.249147045813403, 0.233492536538355,
		0.203167426723066, 0.160078328543346, 0.106939325995318, 0.047175336386512 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder21()
{
	const int n = 11;
	double x[n] = { -0.978228658146057, -0.8870625997680953, -0.7301520055740493, -0.5190961292068118, -0.2695431559523449,
		 0.0000000000000000, 0.2695431559523449, 0.5190961292068118, 0.7301520055740493, 0.8870625997680953, 0.978228658146057 };
	double w[n] = { 0.5566856711617367e-1, 0.1255803694649046, 0.1862902109277342, 0.2331937645919904,
		 0.2628045445102467, 0.2729250867779006, 0.2628045445102467, 0.2331937645919904, 0.1862902109277342, 0.1255803694649046,
		 0.5566856711617367e-1 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder19()
{
	const int n = 10;
	double x[n] = { -0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312,
		 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717 };
	double w[n] = { 0.6667134430868814e-1, 0.1494513491505806, 0.2190863625159820, 0.2692667193099964,
		 0.2955242247147529, 0.2955242247147529, 0.2692667193099964, 0.2190863625159820, 0.1494513491505806,
		 0.6667134430868814e-1 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder15()
{
	const int n = 8;
	double x[n] = { -0.9602898564975362, -0.7966664774136267, -0.5255324099163289, -0.1834346424956498, 0.1834346424956498,
		 0.5255324099163289, 0.7966664774136267, 0.9602898564975362 };
	double w[n] = { 0.1012285362903762, 0.2223810344533745, 0.3137066458778873, 0.3626837833783619,
		 0.3626837833783619, 0.3137066458778873, 0.2223810344533745, 0.1012285362903762 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder9()
{
	const int n = 5;
	double x[n] = { -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640, };
	double w[n] = { 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder7()
{
	const int n = 4;
	double x[n] = { -sqrt((3 + 2 * sqrt(6.0 / 5)) / 7.0), -sqrt((3 - 2 * sqrt(6.0 / 5)) / 7.0), sqrt((3 - 2 * sqrt(6.0 / 5)) / 7.0), sqrt((3 + 2 * sqrt(6.0 / 5)) / 7.0) };
	double w[n] = { (18 - sqrt(30.0)) / 36, (18 + sqrt(30.0)) / 36, (18 + sqrt(30.0)) / 36, (18 - sqrt(30.0)) / 36 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder5()
{
	const int n = 3;
	double x[n] = { -0.7745966692414833, 0.0, 0.7745966692414833 };
	double w[n] = { 0.5555555555555556, 0.8888888888888889, 0.5555555555555556 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline vector<QuadratureNode> SegmentGaussOrder3()
{
	const int n = 2;
	double x[n] = { -0.5773502691896258, 0.5773502691896258 };
	double w[n] = { 1.0, 1.0 };

	vector<QuadratureNode> nodes;

	for (int i = 0; i < n; i++)
		nodes.push_back(QuadratureNode(Point((x[i] + 1) / 2, 0), 0.5 * w[i]));

	return nodes;
}

inline void MultMatrixVector(vector<vector<double>>& A, vector<double>& b, vector<double>& res)
{
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < b.size(); j++)
			res[i] += A[i][j] * b[j];
}

inline void MultMatrixMatrix(vector<vector<double>>& A, vector<vector<double>>& B, vector<vector<double>>& res)
{
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < B.size(); j++)
		{
			res[i][j] = 0;
			for (int k = 0; k < A[i].size(); k++)
				res[i][j] += A[i][k] * B[k][j];
		}
}

inline void MultMatrixMatrix2(vector<vector<double>>& A, vector<vector<double>>& B, vector<vector<double>>& res)
{
	for (int i = 0; i < A.size(); i++)
		for (int j = 0; j < B.size(); j++)
		{
			for (int k = 0; k < B[j].size(); k++)
			{
				res[i][k] = 0;
				res[i][k] += A[i][j] * B[j][k];
			}
		}
}

inline double Gauss18(function<double(double)> f)
{
   double sum = 0;

   for (auto& node : SegmentGaussOrder19())
   {
      double x = node.Node.r;

      sum += node.Weight * f(x);
   }

   return sum;
}

inline double NewtonCotes(double a, double b, function<double(double)> f)
{
   double h = fabs(b - a) / 1000;
   double result = f(a);

   for (double x = h; x < b; x += h)
      result += 2 * f(x);

   result += f(b);
   result *= 7;

   for (double x = a; x < b; x += h)
      result += 32 * f(x + h / 4) + 12 * f(x + h / 2) + 32 * f(x + 3 * h / 4);

   result = result * 0.5 * h / 45;

   return result;
}

inline void Gauss(vector<vector<double>> A, vector<double>& x, vector<double> b)
{
	int N = A.size();
	for (int k = 0; k < N - 1; k++)
	{
		// Обнуление k-ого столбца
		for (int i = k + 1; i < N; i++)
		{
			double t = A[i][k] / A[k][k];
			b[i] -= t * b[k];
			for (int j = k + 1; j < N; j++)
				A[i][j] -= t * A[k][j];
		}
	}

	// Вычисление вектора x
	x[N - 1] = b[N - 1] / A[N - 1][N - 1];
	for (int k = N - 2; k >= 0; k--)
	{
		double sum = 0;
		for (int j = k + 1; j < N; j++)
			sum += A[k][j] * x[j];

		x[k] = (b[k] - sum) / A[k][k];
	}
}