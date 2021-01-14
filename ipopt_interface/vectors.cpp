#include "vectors.h"

vectord operator +(const vectord &v1, const vectord &v2)
{
	vectord out = v1;
	auto e2 = v2.begin();
	for (auto e = out.begin(); e != out.end(); e++, e2++)
		(*e) += (*e2);
	return out;
}

vectord operator -(const vectord &v1, const vectord &v2)
{
	vectord out = v1;
	auto e2 = v2.begin();
	for (auto e = out.begin(); e != out.end(); e++, e2++)
		(*e) -= (*e2);
	return out;
}

vectord operator *(const vectord &v1, const vectord &v2)
{
	vectord out = v1;
	auto e2 = v2.begin();
	for (auto e = out.begin(); e != out.end(); e++, e2++)
		(*e) *= (*e2);
	return out;
}

vectord operator /(const vectord &v1, const vectord &v2)
{
	vectord out = v1;
	auto e2 = v2.begin();
	for (auto e = out.begin(); e != out.end(); e++, e2++)
		(*e) /= (*e2);
	return out;
}

vectord operator /(const vectord &v1, const double &d)
{
	vectord out = v1;
	for (auto e = out.begin(); e != out.end(); e++)
		(*e) /= d;
	return out;
}

vectord operator *(const vectord &v1, const double &d)
{
	vectord out = v1;
	for (auto e = out.begin(); e != out.end(); e++)
		(*e) *= d;
	return out;
}

double  max_abs(const vectord &v1)
{
	double maxx = 0;
	for (auto e = v1.begin(); e != v1.end(); e++)
	{
		if (abs(*e) > maxx)
			maxx = abs(*e);
	}
	return maxx;
}

double interpolate_2d(vectord square_points_vals, double x0, double x1, double y0, double y1, double x, double y)
{
	double xd = (x - x0) / (x1 - x0);
	double yd = (y - y0) / (y1 - y0);
	double c0 = square_points_vals[bin2dec(0, 0)] * (1 - xd) + square_points_vals[bin2dec(1, 0)] * xd;
	double c1 = square_points_vals[bin2dec(0, 1)] * (1 - xd) + square_points_vals[bin2dec(1, 1)] * xd;
	return c0 * (1 - yd) + c1 * yd;
}

double interpolate_3d(vectord cube_points_vals, double x0, double x1, double y0, double y1, double z0, double z1, double x, double y, double z)
{
	double xd = (x - x0) / (x1 - x0);
	double yd = (y - y0) / (y1 - y0);
	double zd = (z - z0) / (z1 - z0);
	double c00 = cube_points_vals[bin2dec(0, 0, 0)] * (1 - xd) + cube_points_vals[bin2dec(1, 0, 0)] * xd;
	double c01 = cube_points_vals[bin2dec(0, 0, 1)] * (1 - xd) + cube_points_vals[bin2dec(1, 0, 1)] * xd;
	double c10 = cube_points_vals[bin2dec(0, 1, 0)] * (1 - xd) + cube_points_vals[bin2dec(1, 1, 0)] * xd;
	double c11 = cube_points_vals[bin2dec(0, 1, 1)] * (1 - xd) + cube_points_vals[bin2dec(1, 1, 1)] * xd;

	double c0 = c00 * (1 - yd) + c10 * yd;
	double c1 = c01 * (1 - yd) + c11 * yd;

	return c0 * (1 - zd) + c1 * zd;

}