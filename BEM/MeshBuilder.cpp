#include "MeshBuilder.h"

void InputGrid(string file, Interval& intervalX, Interval& intervalY)
{
	ifstream in(file);

	in >> intervalX.beg >> intervalX.end >> intervalX.n;
	in >> intervalY.beg >> intervalY.end >> intervalY.n;

	in.close();
}

void BuildRectMesh(Interval& interval, vector<double>& coord)
{
	double beg = interval.beg;
	double end = interval.end;
	int n = interval.n;

	double h = (end - beg) / n;

	for (int i = 0; i <= n; i++)
		coord.push_back(beg + h * i);
}

void BuildCircleMesh(Interval& interval, vector<double>& coord)
{
}

// Создание точек
void CreatePoints(vector<Point>& points, vector<double>& x, vector<double>& y)
{
	int kx = x.size();
	int ky = y.size();

	for (int i = 0; i < kx - 1; i++)
		points.push_back(Point(x[i], y[0]));

	for (int i = 0; i < ky - 1; i++)
		points.push_back(Point(x[kx - 1], y[i])); 

	for (int i = kx - 1; i > 0; i--)
		points.push_back(Point(x[i], y[ky - 1]));

	for (int i = ky - 1; i > 0; i--)
		points.push_back(Point(x[0], y[i]));
}

void CreateRectBEMElements(vector<BEMElement>& elements, int pointSize)
{
	int p = 0;

	for (; p < pointSize - 1; p++)
	{
		BEMElement el;
		el.verts.resize(2);

		el.verts[0] = p;
		el.verts[1] = p + 1;
		el.num = p;

		elements.push_back(el);
	}

	BEMElement el;
	el.verts.resize(2);

	el.verts[0] = p;
	el.verts[1] = elements[0].verts[0];
	el.num = p;

	elements.push_back(el);
}