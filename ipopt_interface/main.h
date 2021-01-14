#pragma once

#include "IpIpoptApplication.hpp"

class Solver;
class MOS_Model_DB;
using namespace std;
//using namespace Ipopt;

//represent parameter binding of two devices
struct Binding
{
	Device *dev1;
	Device *dev2;
	int param1;
	int param2;
	double factor;
};

//An extension for devices to use with GPOptimizer
struct Device_ext
{
	Device *dev;
	int x_index;		//index of first device parameter in optimization parameters array
	vector<int*> px_indices;
	vector<int> x_indices;
	vector<double> x_factors;
	monomial VDS;
	monomial VSB;
	monomial ID;
	monomial mVal;
	vector<monomial*> GPstate;		//device monomial state
};

//Monomial structure for a term in constraint or goal posynomial
//It represent monomial function of device performance parameters which in turn are monomial of device size parameters
struct posterm : monomial
{
	//pointers to devices for each super monomial factors
	vector<Device_ext*> devices;
	//location of first device variable for each super monomial factors
	//foreach device in "devices", the corresponding param_index represent its parameter (gm, gds, vov) used in the factor
	vector<int> param_index;
	//override eval and eval_grad of monomial
	double eval(vectord &x);
	vectord eval_grad(vectord &x);
	vectord eval_hess(vectord &x);
	void Build(smonomial, vector<Device_ext*>*);
private:
	static monomial w;
	static monomial l;
	static monomial i;
	static monomial n;

	int getparamMap(string mod, string pname);
	monomial* getMonomial_Model(Device_ext *device, int factor_index);
};

//Multi terms of posterm(s) added together
struct posynomial
{
	vector<posterm> terms;
	//temp store the final value of goals and constraints after finishing ipopt optimization
	double final_value;			//It stores the last evaluated value
	string name;
	double overdesign_factor = 1;
	double eval(vectord &x);
	vectord eval_grad(vectord &x);
	vectord eval_hess(vectord &x, double multiplier);
	vectord eval_hess(vectord &x, const vectord &add, double multiplier);
	//void eval_grad(vectord &x, vectord &out);
	//void eval_hess(vectord &x, const vectord &add, double multiplier, vectord &out);
};

//abstract optimizer
class Optimizer
{
protected:
	vector<Binding> bindings;
public:
	Circuit * CDB;
	virtual void add_constraint(vectorst components, int type, string *name = NULL) = 0;
	virtual void set_goal(vectorst components) = 0;
	virtual int map_param(string pName);
	virtual void add_binding(string name1, string name2, string factor = "");
	virtual void run(int maxIter) = 0;
	virtual void run_test(int maxIter) = 0;
	inline Optimizer(Circuit *cdb) : CDB{ cdb } {};
	virtual ~Optimizer() {};
};

//Geometric programming optimizer
class GPOptimizer : public Optimizer
{
	vector<posynomial> constraints;
	vector<int> constraints_types;
	posynomial goal;
	MOS_Model_DB *mmdb = NULL;
	Solver *slv = NULL;
	double update_factor;			//solver update factor
	vector<Device_ext*> dev_ext;		//Constructed upon object creation except for their GPstate
	void initialize_app(Ipopt::SmartPtr<Ipopt::IpoptApplication> &app);
	Ipopt::SmartPtr<Ipopt::TNLP> ipm;
	vector<range> initialize_devices();
	void wite_solution_back(vectord& solution);
	void update_gp_models(vectord& indeps, bool initialize);
	bool evaluate_corners(vector<double> &typ_res);
	double compare_spice(vectord solution);
public:
	void add_constraint(vectorst posynomial, int type, string *name = NULL);
	void set_goal(vectorst posynomial);
	void run(int maxIter);
	void specs_correction(vectord &indeps, vectord &sim_vals);
	void run_test(int maxIter);
	void BuildModels(string filepath, vectorst ranges);
	void LoadModels(string filepath);
	GPOptimizer(Circuit *cdb);
	~GPOptimizer();
};
