#pragma once

struct Point
{
	Point(double r, double z) : r(r), z(z) {}
	Point() : r(0.0), z(0.0) {}
	double r, z;

	Point operator+(Point a)
	{
		return Point(r + a.r, z + a.z);
	}

	Point operator-(Point a)
	{
		return Point(r - a.r, z - a.z);
	}

	double operator*(Point a)
	{
		return r * a.r + z * a.z;
	}

	Point operator*(double constant)
	{
		return Point(r * constant, z * constant);
	}

	Point operator/(double constant)
	{
		return Point(r / constant, z / constant);
	}

	bool operator==(Point a)
	{
		return (r == a.r && z == a.z);
	}

	double Norm()
	{
		return sqrt(Point(r, z) * Point(r, z));
	}

	Point Normalize()
	{
		Point v(r, z);
		return v / v.Norm();
	}
};