#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED
#include <string>
#include <vector>
#include <math.h>

using namespace std;

//General independent type definitions

//general parameter for .model, .param, instants
struct parameter
{
	const string name;	//name is set one time only when the struct is constructed
	double value;
	static vector<parameter*>* allParams;

	//should convert eng notation to sci
	//should split expression first 1u+vgs-y ==> 1 u + vgs - y
	parameter(string p_name, bool save);
	parameter(string p_name, string &expression, bool save);
	parameter(string p_name, double p_value, bool save);

	static double try_parsing(string expression);

	static double get_value_of(string p_name);

	static parameter* isParam(string p_name)
	{
		for (auto param = (*allParams).begin(); param != (*allParams).end(); param++)
			if ((**param).name == p_name)
				return &**param;
		return NULL;
	}

	static string sci2eng(double number);

	void old_eng2sci(string &expression)
	{
		char unit = *(expression.end() - 1);
		if (unit >= '0' && unit <= '9')
			return;		//last char is numerical
		expression.erase(expression.end() - 1);
		string sci;
		switch (unit)
		{
		case 'k':
			sci = "e+3"; break;
		case 'm':
			sci = "e-3"; break;
		case 'u':
			sci = "e-6"; break;
		case 'n':
			sci = "e-9"; break;
		case 'p':
			sci = "e-12"; break;
		}
		expression.insert(expression.end(), sci.begin(), sci.end());
	}
	//Remove the copy constructor to prevent copying
	parameter operator =(const parameter &) { return parameter(name, value, false); };
private:
	static bool initialized;
	static void initialize();
};

typedef class vector<string> vectorst;

//typedef class Ngsimulator_dl Ngsimulator;

#endif // DEFS_H_INCLUDED

/*
Need to create a multidimensional array template
it is based on one-dimension array with indexes
1st-dim-index + (#elem_in_1st-dim)* 2nd-dim-index + (#elem_in_1st-dim)*(#elm_in_2nd-dim) * 3rd-dim-index + ....
#elements of each dimension should be fixed and can be different

*/
