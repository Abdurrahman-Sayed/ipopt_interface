#pragma once
#include <string>
#include <vector>
#include <math.h>
using namespace std;

typedef class vector<double> vectord;

//vector of doubles operations
vectord operator +(const vectord &v1, const vectord &v2);
vectord operator -(const vectord &v1, const vectord &v2);
vectord operator *(const vectord &v1, const vectord &v2);
vectord operator /(const vectord &v1, const vectord &v2);
vectord operator /(const vectord &v1, const double &d);
vectord operator *(const vectord &v1, const double &d);
double  max_abs(const vectord &v1);
inline const int bin2dec(const int b0, const int b1, const int b2)
{
	return b0 + 2 * b1 + 4 * b2;
}
inline const int bin2dec(const int b0, const int b1)
{
	return b0 + 2 * b1;
}
double interpolate_2d(vectord square_points_vals, double x0, double x1, double y0, double y1, double x, double y);
//cube_points_vals are sorted as b0=x , b1=y , b2=z
double interpolate_3d(vectord cube_points_vals, double x0, double x1, double y0, double y1, double z0, double z1 , double x, double y, double z);