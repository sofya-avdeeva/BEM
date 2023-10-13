#pragma once

#include <vector>
#include <string>
#include <fstream>

#include "Point.h"

using namespace std;

struct Interval
{
	double beg, end;
	int n;
};

struct BEMElement
{
	vector<int> verts;
	int num;

	bool operator==(BEMElement e)
	{
		return (verts[0] == e.verts[0] && verts[1] == e.verts[1]);
	}

	bool operator!=(BEMElement e)
	{
		return (verts[0] != e.verts[0] && verts[1] != e.verts[1]);
	}
};

void InputGrid(string file, Interval& intervalX, Interval& intervalY);
void BuildRectMesh(Interval& interval, vector<double>& coord);
void BuildCircleMesh(Interval& interval, vector<double>& coord);
void CreatePoints(vector<Point>& points, vector<double>& x, vector<double>& y);
void CreateRectBEMElements(vector<BEMElement>& elements, int pointSize);