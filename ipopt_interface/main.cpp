#include <iostream>
#include <fstream>
#include "main.h"
#include "ipopt.h"
#include "optimizergp.h"

#include <chrono>
std::ofstream Log;
using namespace std;

monomial posterm::w = monomial(1, { 1, 0 , 0 });
monomial posterm::l = monomial(1, { 0, 1 , 0 });
monomial posterm::i = monomial(1, { 0, 0 , 1 });
monomial posterm::n = monomial(1, { 1 });
double shrink_factor = 100;

chrono::microseconds eval_clock;
chrono::microseconds grad_clock;
chrono::microseconds hess_clock;

///#define MEAS_TIME
///#define DEBUG_DEVICES

//Future work:
//.match m1[id] m5[id] 0.5		:    m1[id] = m5[id] * 0.5	, m5[id] is the reference

int posterm::getparamMap(string modType, string pname)
{
	if (modType == "nmos" || modType == "pmos")
	{
		if (pname == "vov")		//It is allowed to use veff (vov) as part of optimization goal/constraint
			return 0;
		if (pname == "gds")
			return 1;
		if (pname == "gm")
			return 2;
		if (pname == "vth")
			return 3;
		if (pname == "cgg")
			return 4;
		if (pname == "cdd")
			return 5;
		if (pname == "css")
			return 6;
		if (pname == "gmb")
			return 7;
		if (pname == "w")
			return 8;
		if (pname == "l")
			return 9;
		if (pname == "vds")
			return 10;
		if (pname == "vsb")
			return 11;
		if (pname == "id")
			return 12;
		//Need to add more parameters from const (id, vds). may use negative values for indices of const parameters
	}
	return -10000;	//return invalid output for invalid input
}

monomial* posterm::getMonomial_Model(Device_ext *device, int factor_index)
{
	int pi = param_index[factor_index];
	if (pi < 0)
		return &device->mVal;
	else if (pi <= 7)
		return (device->GPstate)[pi];
	else if (pi == 8)
		return &w;
	else if (pi == 9)
		return &l;
	else if (pi == 10)
		return &device->VDS;
	else if (pi == 11)
		return &device->VSB;
	else if (pi == 12)
		//	return &device->ID;
		return &i;
	else
		return NULL;
}

void GPOptimizer::add_constraint(vectorst stposynomial, int type, string *name)
{
	if (type == 0 && stposynomial.size() > 1)
		throw new exception("equality constraints must have one term (monomial)");
	posynomial p;
	for (int i = 0; i < stposynomial.size(); i++)
	{
		smonomial m(stposynomial[i]);
		posterm pm;
		pm.Build(m, &dev_ext);
		p.terms.push_back(pm);
	}
	if (name != NULL)
		p.name = *name;
	constraints.push_back(p);
	constraints_types.push_back(type);
}

void GPOptimizer::set_goal(vectorst stposynomial)
{
	for (int i = 0; i < stposynomial.size(); i++)
	{
		smonomial m(stposynomial[i]);
		posterm pm;
		pm.Build(m, &dev_ext);
		goal.terms.push_back(pm);
	}
}

void posterm::Build(smonomial m, vector<Device_ext*> *devs)
{
	coef = m.coef;
	powers = move(m.powers);
	for (auto n : m.names)
	{
		size_t br = n.find_first_of('[');
		string devName;
		string paramName;
		if (br == string::npos)
		{
			devName = n;
			paramName = "";
		}
		else
		{
			devName = n.substr(0, br);
			paramName = n.substr(br + 1, n.size() - br - 2);
		}
		for (auto d = devs->begin(); d != devs->end(); ++d)
		{
			if ((*d)->dev->Name == devName)
			{
				devices.push_back((*d));
				if ((*d)->dev->model != NULL)
					param_index.push_back(getparamMap((*d)->dev->model->type, paramName));
				else
					param_index.push_back(-1);
				break;
			}
		}
	}
}

void GPOptimizer::BuildModels(string filepath, vectorst ranges)
{
	//Build all models databases
	//MOS models

	mmdb->Build_DB(filepath, ranges);
	//mmdb->Build_DB(rangeN(0.0001, 0.001, 20), rangeN(1.5, 3.3, 19), rangeN(0.0, 0.5, 6), rangeN(10e-6, 100e-6, 6), rangeN(0.6e-6, 10e-6, 5), filepath);
}

void GPOptimizer::LoadModels(string filepath)
{
	//Load models databases from file
	mmdb->Load_DB(filepath);
}

int Optimizer::map_param(string pName)
{
	if (pName == "w")
		return 0;
	if (pName == "l")
		return 1;
	if (pName == "id")
		return 2;
	return -1;
}

void Optimizer::add_binding(string name1, string name2, string factor)
{
	Binding b;
	double fac = 1;
	if (factor != "")
		fac = stod(factor);
	b.factor = fac;

	string devName1, devName2;
	size_t br = name1.find_first_of('[');
	if (br == string::npos)
	{
		devName1 = name1;
		b.param1 = -1;
	}
	else
	{
		devName1 = name1.substr(0, br);
		string paramName = name1.substr(br + 1, name1.size() - br - 2);
		b.param1 = map_param(paramName);
	}

	br = name2.find_first_of('[');
	if (br == string::npos)
	{
		devName2 = name2;
		b.param2 = -1;
	}
	else
	{
		devName2 = name2.substr(0, br);
		string paramName = name2.substr(br + 1, name2.size() - br - 2);
		b.param2 = map_param(paramName);
	}

	for (auto d = CDB->Devices.begin(); d != CDB->Devices.end(); ++d)
	{
		if ((*d)->Name == devName1)
			b.dev1 = *d;
		if ((*d)->Name == devName2)
			b.dev2 = *d;
	}
	bindings.push_back(b);
}

GPOptimizer::GPOptimizer(Circuit *cdb) : Optimizer(cdb)
{
	//mmdb = new MOS_Model_DB(20, 10, 5, 10, 5);
	//mmdb = new MOS_Model_DB(3, 2, 2, 3, 2);
	//mmdb = new MOS_Model_DB(1, 1, 1, 1, 1);
	mmdb = new MOS_Model_DB_T();
	//Create GPmodels based on SPICE models (upgrade)
	for (auto m = cdb->Models.begin(); m != cdb->Models.end(); m++)
	{
		GPmodel *gpm = NULL;
		if ((*m)->type == "nmos" || (*m)->type == "pmos")
		{
			gpm = mmdb->createGPmodel(dynamic_cast<spice_model*>(*m));
			*m = gpm;
		}
		//Regard only upgradable Models
		if (gpm == NULL)
			continue;
		//Correct model pointers in devices in circuit database to point to the upgraded version GPmodel
		for (auto d : cdb->Devices)
		{
			if (!d->model)
				continue;
			if (!d->optimize)
				continue;
			if (d->model->Name == (*m)->Name)
			{
				//The current device has the current model
				Device_ext *de = new Device_ext();
				de->dev = d;
				de->VDS = monomial(1, {});
				de->VSB = monomial(1, {});
				de->ID = monomial(1, {});
				size_t s = gpm->names.size();
				de->x_indices.resize(s);
				de->px_indices.resize(s);
				de->x_factors.resize(s);
				for (int i = 0; i < s; i++)
				{
					de->px_indices[i] = &de->x_indices[i];
					de->x_factors[i] = 1;
				}
				dev_ext.push_back(de);
				d->model = gpm;
			}
		}
	}
	//Add devices that don't have models
	for (auto d : cdb->Devices)
	{
		if (d->model)
			continue;
		if (!d->optimize)
			continue;
		Device_ext *de = new Device_ext();
		de->dev = d;
		de->mVal = monomial(1, { 1 });
		de->x_indices.resize(1);
		de->px_indices.resize(1);
		de->x_factors.resize(1);
		de->px_indices[0] = &de->x_indices[0];
		de->x_factors[0] = 1;
		dev_ext.push_back(de);
	}
	//add sources that need to be optimized
	for (auto d : cdb->Sources)
	{
		if (!d->optimize)
			continue;
		Device_ext *de = new Device_ext();
		de->dev = d;
		de->mVal = monomial(1, { 1 });
		de->x_indices.resize(1);
		de->px_indices.resize(1);
		de->x_factors.resize(1);
		de->px_indices[0] = &de->x_indices[0];
		de->x_factors[0] = 1;
		dev_ext.push_back(de);
	}
}

GPOptimizer::~GPOptimizer()
{
	delete mmdb;
	for (auto f : dev_ext)
		delete f;
}

vector<range> GPOptimizer::initialize_devices()
{
	vector<range> ranges;
	size_t c_index = 0;
	for (auto de = dev_ext.begin(); de != dev_ext.end(); de++)
	{
		GPmodel *gpm = dynamic_cast<GPmodel*>((*de)->dev->model);

		//Handle bindings
		Device_ext* ref_dev = NULL;
		Binding *mb;
		vector<int> refs;
		for (auto b = bindings.begin(); b != bindings.end(); b++)
		{
			if ((*de)->dev == (*b).dev1)
			{
				mb = &*b;
				//search for its binded reference
				for (auto rd = dev_ext.begin(); rd != dev_ext.end(); rd++)
					if ((*rd)->dev == (*b).dev2)
						ref_dev = *rd;
				//detected match dev[j] dev[ref_index]



				if (mb->param1 == -1)
				{
					for (int i = 0; i < gpm->names.size(); i++)
					{
						(*de)->px_indices[i] = &ref_dev->x_indices[i];
						(*de)->x_factors[i] = mb->factor;
						refs.push_back(i);
					}
				}
				else
				{
					(*de)->px_indices[mb->param1] = &ref_dev->x_indices[mb->param2];
					(*de)->x_factors[mb->param1] = mb->factor;
					refs.push_back(mb->param1);
				}
			}
		}
		if (ref_dev == NULL)
		{
			(*de)->x_index = (int)c_index;	//non-referenced device
			for (int i = 0; i < (*de)->x_indices.size(); i++)
			{
				(*de)->x_indices[i] = (int)c_index++;
			}
			//if (gpm == NULL)
			//	c_index += 1;
			//else
			//	c_index += gpm->names.size();
		}
		else
		{
			//de->x_index = -(int)c_index;	//referenced device
			//dev_ext[device_index]
			//	---> x_index is first index of first indep parameter (w, l) of cdb->Devices[device_index] 
			//        in all indep parameters vector
			//	---> GPstate is state data of device cdb->Devices[device_index]
			vector<int> non_refs;
			for (int i = 0; i < (*de)->x_indices.size(); i++)
			{
				bool found = false;
				for (int j = 0; j < refs.size(); j++)
					if (refs[j] == i)
						found = true;
				if (!found)
					non_refs.push_back(i);
			}
			for (int i = 0; i < non_refs.size(); i++)
			{
				(*de)->x_indices[non_refs[i]] = (int)c_index++;
			}
			//continue;
		}
		//Initialize design parameter ranges for each device (only for unbinded devices)
		if (gpm)
		{
			if (gpm->type == "nmos" || gpm->type == "pmos")
			{
				MOS *d = dynamic_cast<MOS*>((*de)->dev);
				//add ranges only to non-referenced parameters
				int refsh1 = 0;
				for (int k = 0; k < refs.size(); k++)
					refsh1 += 1 << refs[k];
				if (!(refsh1 & 1))
					ranges.push_back(range(log(gpm->independent_ranges[0].begin), log(gpm->independent_ranges[0].end)));
				if (!(refsh1 & 2))
					//ranges.push_back(range(log(2e-6), log(2e-6)));
					ranges.push_back(range(log(gpm->independent_ranges[1].begin), log(gpm->independent_ranges[1].end)));
				if (!(refsh1 & 4))
					ranges.push_back(range(log(1e-7), log(1e-2)));
			}
		}
		else
		{
			if (refs.size() > 0)
				continue;
			//switch ideal device type
			switch ((*de)->dev->Name[0])
			{
			case 'r':
				ranges.push_back(range(log(1e1), log(1e7)));
				break;
			case 'c':
				ranges.push_back(range(log(1e-14), log(10e-9)));		//cap range. Need to be revised
				break;
			case 'v':
				ranges.push_back(range(log(0.01), log(5)));
				break;
			case 'i':
				ranges.push_back(range(log(1e-7), log(1e-2)));
				break;
			}
		}
	}
	if (ranges.size() != c_index)
		throw new exception("Inconsistence devices and/or model");
	return ranges;
}

void GPOptimizer::initialize_app(Ipopt::SmartPtr<Ipopt::IpoptApplication> &app)
{
	app->Options()->SetNumericValue("tol", 1e-5);
	//app->Options()->SetNumericValue("constr_viol_tol", 1);
	//app->Options()->SetNumericValue("dual_inf_tol", 2);
	app->Options()->SetStringValue("recalc_y", "yes");
	//app->Options()->SetNumericValue("recalc_y_feas_tol", 0.001);
	//app->Options()->SetStringValue("alpha_for_y", "min-dual-infeas");
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("bound_mult_init_method", "mu-based");
	//app->Options()->SetStringValue("mehrotra_algorithm", "yes");
	app->Options()->SetIntegerValue("max_iter", 15);
	app->Options()->SetNumericValue("mu_init", 0.5);
	app->Options()->SetIntegerValue("print_level", 5);

	//app->Options()->SetNumericValue("derivative_test_perturbation", 1e-9);
	//app->Options()->SetStringValue("derivative_test", "first-order");
	//app->Options()->SetStringValue("derivative_test", "second-order");
}

void GPOptimizer::run(int maxIter)
{
	if (slv == NULL)
	{
		slv = new Simulator(CDB);
		////Step 1 : Initialize all variables
		//initialize circuit solver
		slv->Initialize();
	}
	//initialize devices' independant parameters index mapping and ranges
	vector<range> ranges = initialize_devices();
	size_t numvars = ranges.size();
	vectord indeps(numvars);		//give initial point to IPM optimizer and take final point
	vectord old_indeps(numvars);		//store old point to		
	vector<double> sim_vals;		//auxilary variable to take simulated specs from the simulator
									///////Step 2: Update GP models of all devices
	update_gp_models(indeps, true);		//this is only to collect info about scaling... Not anymore
										//initialize IPM solver
	ipm = new IPMOpt(&indeps, &ranges, &goal, &constraints, &constraints_types);
	/*
	IPMOpt *gpipm = dynamic_cast<IPMOpt*>(GetRawPtr<TNLP>(ipm));
	//set a custom stop function for the IPMOpt
	gpipm->custom_stop_func = &GPOptimizer::change_stop;

	bool (GPOptimizer::*custom_stop_func)(const vectord& initial) = NULL;
	*/
	update_factor = 0.3;

	Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
	initialize_app(app);
	app->Options()->SetNumericValue("obj_scaling_factor", 1 / goal.eval(indeps));
	// Initialize the IpoptApplication and process the options
	Ipopt::ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Ipopt::ApplicationReturnStatus::Solve_Succeeded)
	{
		std::cout << std::endl << std::endl << "*** Error during initialization of IPM solver!" << std::endl;
		return;
	}
	int state = 10;	//convergence indicator  0:complete conv    1:close to conv   10:far from conv    -1:simulation didn't converge
	int old_state = state;
	bool corners_optim_phase = false;
#ifdef MEAS_TIME
	auto start = chrono::high_resolution_clock::now();
	auto total_start = chrono::high_resolution_clock::now();
	auto stop = chrono::high_resolution_clock::now();
	chrono::microseconds gp_clock(chrono::duration_values<long long>::zero());
	chrono::microseconds sim_clock(chrono::duration_values<long long>::zero());
	chrono::microseconds get_model_clock(chrono::duration_values<long long>::zero());
	chrono::microseconds compare_clock(chrono::duration_values<long long>::zero());
	chrono::microseconds total_time(chrono::duration_values<long long>::zero());
#endif
	//Main optimization loop
	for (int i = 1; i <= maxIter; i++)
	{
		//app->Options()->SetIntegerValue("max_iter", 110 - state*10);		//from 2 to 11
		old_indeps = indeps;
		//update_gp_models(indeps, false);
		///////Step 3: Solve GP
		//app->Options()->SetIntegerValue("max_iter", state * 10 + 10);

#ifdef MEAS_TIME
		start = chrono::high_resolution_clock::now();
#endif

		if (i == 1 || true)
			status = app->OptimizeTNLP(ipm);
		else
			status = app->ReOptimizeTNLP(ipm);


#ifdef MEAS_TIME
		stop = chrono::high_resolution_clock::now();
		gp_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif

		int ckt_converged = 0;
		while (ckt_converged <= 0)
		{

#ifdef MEAS_TIME
			start = chrono::high_resolution_clock::now();
#endif

			//Write the solution to devices
			wite_solution_back(indeps);
			old_state = state;

			///////Step 4: Solve the circuit with the new sizing (Simulation Loop)
			state = slv->Update(ckt_converged == -1, update_factor, sim_vals, { "GBW", "phase margin", "Av", "power" });
#ifdef MEAS_TIME
			stop = chrono::high_resolution_clock::now();
			sim_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
			start = chrono::high_resolution_clock::now();
#endif

			///////Repeat Step 2
			update_gp_models(indeps, false);
#ifdef MEAS_TIME
			stop = chrono::high_resolution_clock::now();
			get_model_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif
			//Correct GBW and phase margin
			specs_correction(indeps, sim_vals);

			if (state != -1)
			{
				if (ckt_converged == -1)
					state = 10; //converged but not in the first trial
				ckt_converged = 1;
			}
			else
			{
				indeps = (indeps + old_indeps) / 2;
				ckt_converged = -1;
			}
		}

		//double factor = 0.2 / 9 * state + 0.9;
		//update_factor /= factor;
		//shrink_factor *= factor;
		//if (update_factor > 1.0)
		//	update_factor = 1.0;

		if (status == Ipopt::ApplicationReturnStatus::Solve_Succeeded)
		{
			//if (shrink_factor > 5)
			//	shrink_factor = 5;
			//shrink_factor *= 0.5 * state / 10;
			//update_factor = 1.0;
			update_factor = update_factor;
		}
		else if (status == Ipopt::ApplicationReturnStatus::User_Requested_Stop)
			;//	shrink_factor /= 0.9;
		else if (status == Ipopt::ApplicationReturnStatus::Infeasible_Problem_Detected)
			shrink_factor = shrink_factor;//indeps = indeps * 0.0 + old_indeps * 1.0;		//accept only a small portion of variable change

#ifdef MEAS_TIME
		start = chrono::high_resolution_clock::now();
#endif
		double comp = compare_spice(indeps);
#ifdef MEAS_TIME
		stop = chrono::high_resolution_clock::now();
		compare_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif

		if (comp == 0)
		{
			comp = comp;
			//update_factor = 1;
		}

		shrink_factor = 2 * comp + 0.01;
		std::cout << ((Simulator*)slv)->ac_test();
		if (comp == 0 && state == 0 && i > 3)
		{


#ifdef MEAS_TIME
			stop = chrono::high_resolution_clock::now();
			total_time += chrono::duration_cast<chrono::microseconds>(stop - total_start);

			std::cout << "total_time = " << total_time.count() << endl;
			std::cout << "gp_clock = " << gp_clock.count() << endl;
			std::cout << "sim_clock = " << sim_clock.count() << endl;
			std::cout << "get_model_clock = " << get_model_clock.count() << endl;
			std::cout << "compare_clock = " << compare_clock.count() << endl;
			std::cout << "eval_clock = " << eval_clock.count() << endl;
			std::cout << "grad_clock = " << grad_clock.count() << endl;
			std::cout << "hess_clock = " << hess_clock.count() << endl;
#endif

			if (corners_optim_phase)
				;// break;
				 //typical optimization succeded. optimize for corners
			if (evaluate_corners(sim_vals))
				break;

			if (compare_spice(indeps) == 0)
				break;
			corners_optim_phase = true;
			shrink_factor = 1;
			//update_factor = update_factor / 3;
			//app->Options()->SetNumericValue("mu_init", 0.03);
		}
	}
	std::cout << ((Simulator*)slv)->ac_test();
}

void GPOptimizer::specs_correction(vectord &indeps, vectord &sim_vals)
{
	for (Index i = 0; i < constraints.size(); i++)
	{
		if (constraints[i].name == ":gbw")
		{
			double GBWspec = parameter::get_value_of("gbwspec");
			double c = constraints[i].eval(indeps) / constraints[i].overdesign_factor;
			double correction = GBWspec / (c * sim_vals[0]) * 1.02;
			if (isfinite<double>(correction))
			{
				for (int k = 0; k < constraints[i].terms.size(); k++)
					constraints[i].terms[k].coef *= correction;
			}
		}
		else if (constraints[i].name == ":phase_m")
		{
			double cpm = parameter::get_value_of("cpm2");
			double c = constraints[i].eval(indeps) / constraints[i].overdesign_factor;
			double correction = 0;
			if (sim_vals[1] < 0)
				correction = 2;
			else
				correction = cpm / (c * 1 / tan(3.14159 / 180 * (90 - sim_vals[1])));
			if (correction < 0) correction = 0.1;
			if (isfinite<double>(correction))
			{
				for (int k = 0; k < constraints[i].terms.size(); k++)
					constraints[i].terms[k].coef *= correction;
			}
		}
		else if (constraints[i].name == ":av")
		{
			double cav = parameter::get_value_of("cav");
			double c = constraints[i].eval(indeps) / constraints[i].overdesign_factor;;
			double correction = cav / (c *sim_vals[2]);
			if (isfinite<double>(correction))
			{
				for (int k = 0; k < constraints[i].terms.size(); k++)
					constraints[i].terms[k].coef *= correction;
			}
		}
		else if (constraints[i].name == ":power")
		{
			double cp = parameter::get_value_of("cp");
			double cSpec = parameter::get_value_of("vdd") / cp;
			double c = constraints[i].eval(indeps) / constraints[i].overdesign_factor;;
			double correction = sim_vals[3] / (cSpec * c);
			if (isfinite<double>(correction))
			{
				for (int k = 0; k < constraints[i].terms.size(); k++)
					constraints[i].terms[k].coef *= correction;
			}
		}
	}
}

double GPOptimizer::compare_spice(vectord solution)
{
#ifdef DEBUG_DEVICES
	for (auto r : dev_ext)
	{
		if (r->dev->model)
		{
			GPmodel *m = dynamic_cast<GPmodel*>(r->dev->model);
			if (m == NULL)
				throw new exception("Can't retrieve GPmodel of a device");
			if (r->dev->model->type == "nmos" || r->dev->model->type == "pmos")
			{
				MOS *d = dynamic_cast<MOS*>(r->dev);
				//Set device state (based on const parameters). State is monomials for required parameters calculation
				vector<double> sol = { solution[*(r->px_indices[0])], solution[*(r->px_indices[1])] ,solution[*(r->px_indices[2])] };
				double vov = r->GPstate[0]->leval(sol);
				double gds = r->GPstate[1]->leval(sol) * 1e6;
				double gm = r->GPstate[2]->leval(sol) * 1e6;
				gds = gds;
				if (d->Name == "m2")
					cout << " gds2=" + to_string(gds) + " gm2=" + to_string(gm);
				if (d->Name == "m6")
					cout << " gds6=" + to_string(gds);
				if (d->Name == "m4")
					cout << " gds4=" + to_string(gds);
				if (d->Name == "m7")
					cout << " gds7=" + to_string(gds) + " gm7=" + to_string(gm);
			}
		}
	}
#endif
	cout << "\n";
	//print solution after rounding
	std::cout << std::endl << "After rounding:" << std::endl << "Objective value:  ";
	double goal_val = goal.eval(solution);
	std::cout << "f(x*) = " << goal_val << "   ";
	double cav = parameter::get_value_of("cav");
	std::cout << " goal predicted:  " << 20 * log10(cav / goal_val) << endl;

	std::cout << std::endl << "Final value of the constraints:" << std::endl;
	double ret = 0;
	for (Index i = 0; i < constraints.size(); i++)
	{
		double c = constraints[i].eval(solution);
		std::cout << constraints[i].name << " = " << c;
		if (constraints[i].name.size() > 1 && constraints[i].name[1] == 'i')
		{
			std::cout << std::endl;
			continue;
		}
		if (constraints_types[i] == -1)
			if (c > 1.01)
				ret += c - 1;
		if (constraints_types[i] == 0)
			if (c > 1.1 || c < 0.9)
				ret += abs(c - 1);
		double spec = 1;
		if (constraints[i].name == ":gbw")
		{
			spec = parameter::get_value_of("gbwspec");
			std::cout << "  predicted:  " << constraints[i].overdesign_factor * spec / c * 1e-6 << " MHz";
		}
		//else if (constraints[i].name == ":phase_m")
		//{
		//	spec = parameter::get_value_of("cpm");
		//	std::cout << "  predicted:  " << 90 - atan(spec / c ) * 180 / 3.14159;
		//}
		else if (constraints[i].name == ":sr")
		{
			spec = parameter::get_value_of("sr");
			std::cout << "  predicted:  " << constraints[i].overdesign_factor * spec / c * 1e-6 << "V/us";
		}

		std::cout << std::endl;
	}
	return ret;
}

void GPOptimizer::update_gp_models(vectord& indeps, bool initialize)
{
	int i = 0;
	for (auto r : dev_ext)
	{
		if (r->dev->model)
		{
			GPmodel *m = dynamic_cast<GPmodel*>(r->dev->model);
			if (m == NULL)
				throw new exception("Can't retrieve GPmodel of a device");
			if (r->dev->model->type == "nmos" || r->dev->model->type == "pmos")
			{
				MOS *d = dynamic_cast<MOS*>(r->dev);
				//initialize independant variables from all devices except binded devices
				if (initialize)
				{
					//initialize optimization with device parameters
					indeps[*(r->px_indices[0])] = log(d->W);
					indeps[*(r->px_indices[1])] = log(d->L);
					indeps[*(r->px_indices[2])] = log(d->Current);
				}
				//Set device state (based on const parameters). State is monomials for required parameters calculation
				//r->GPstate = mmdb->get_model(m, d->Current, d->Vds, d->Vsb);
				MOS_Model_DB_T *mmdb_t = dynamic_cast<MOS_Model_DB_T*>(mmdb);
				vectord sim_s;
				vectord *psim_s = NULL;
				if (mmdb_t)
				{
					Simulator* sim = ((Simulator*)slv);
					if (sim && sim->ready() && update_factor == 1.0)
					{
						sim_s = ((Simulator*)slv)->get_device_sim_state(d, mmdb_t->required_params);
						psim_s = &sim_s;
					}
					range wr, lr, ir;
					wr.begin = exp(indeps[*(r->px_indices[0])]) / (MAX_CHANGE_RATIO);
					wr.end = exp(indeps[*(r->px_indices[0])]) * MAX_CHANGE_RATIO;
					lr.begin = exp(indeps[*(r->px_indices[1])]) / (MAX_CHANGE_RATIO);
					lr.end = exp(indeps[*(r->px_indices[1])]) * MAX_CHANGE_RATIO;
					ir.begin = exp(indeps[*(r->px_indices[2])]) / (MAX_CHANGE_RATIO);
					ir.end = exp(indeps[*(r->px_indices[2])]) * MAX_CHANGE_RATIO;
#ifdef INTERPOLATION
					r->GPstate = mmdb_t->get_model_interp(m, d->Current, d->Vds, d->Vsb, wr, lr, ir, i, psim_s);
#else
					r->GPstate = mmdb_t->get_model(m, d->Current, d->Vds, d->Vsb, wr, lr, ir, i, psim_s);
#endif
				}
				else
#ifdef INTERPOLATION
					r->GPstate = mmdb->get_model_interp(m, d->Current, d->Vds, d->Vsb, i, psim_s);
#else
					r->GPstate = mmdb->get_model(m, d->Current, d->Vds, d->Vsb, psim_s);
#endif
				r->VDS.coef = abs(d->Vds);
				r->VSB.coef = abs(d->Vsb);
				r->ID.coef = abs(d->Current);
				//for (int h = 0; h < m->monmodels[0].powers.size(); h++)
			}
			i++;
		}
		else
		{
			if (initialize)
			{
				indeps[*(r->px_indices[0])] = log(r->dev->Value);
				r->mVal.coef = 1;
			}
		}
	}
}
void GPOptimizer::wite_solution_back(vectord& solution)
{
	//Get grid size
	double l = parameter::get_value_of("lambda");
	for (auto r : dev_ext)
	{
		if (r->dev->model)
		{
			if (r->dev->model->type == "nmos" || r->dev->model->type == "pmos")
			{
				MOS *d = dynamic_cast<MOS*>(r->dev);
				d->W = exp(solution[*(r->px_indices[0])]);
				d->L = exp(solution[*(r->px_indices[1])]);
				//for debug only
				std::cout << d->Name << "  W= " << d->W << "  L= " << d->L << "  Id= " << d->Current;
				if (std::isfinite<double>(l))
				{
					d->W = round(d->W / l)*l;
					d->L = round(d->L / l)*l;
					solution[*(r->px_indices[0])] = log(d->W);
					solution[*(r->px_indices[1])] = log(d->L);
					//					std::cout << "    After rounding:  W= " << d->W << "  L= " << d->L;
				}
				std::cout << endl;
				vectord x = { d->W , d->L };
				d->Vov = r->GPstate[0]->leval(x);
				d->Vt = r->GPstate[3]->leval(x);
			}
		}
		else
		{
			r->dev->Value = exp(solution[*(r->px_indices[0])]);
			std::cout << r->dev->Name << "  = " << r->dev->Value << endl;
		}
	}
	//print solution after rounding
	//	std::cout << std::endl << "After rounding:" << std::endl << "Objective value:  ";
	//	std::cout << "f(x*) = " << goal.final_value << std::endl;
	//	std::cout << std::endl << "Final value of the constraints:" << std::endl;
	//	for (Index i = 0; i < constraints.size(); i++)
	//	{
	//		std::cout << "g(" << i << ") = " << constraints[i].eval(solution) << std::endl;
	//	}
}

bool GPOptimizer::evaluate_corners(vector<double> &typ_res)
{
	//save final value of typical corner
	//vector<double> typical;
	//for (Index i = 0; i < constraints.size(); i++)
	//	typical.push_back(constraints[i].final_value);
	//	vectorst const_names = { ":gbw" , ":phase_m", ":av", ":power" };
	bool all_satisfied = true;
	vectord cor_fac(typ_res.size(), 0.5);	//factor less than 1 is not allowed
	vectord worst_res(typ_res.size());
	//save old overdesign factors
	for (Index i = 0; i < constraints.size(); i++)
		constraints[i].final_value = constraints[i].overdesign_factor;
	for (auto cor = CDB->sCorners.begin(); cor != CDB->sCorners.end(); cor++)
	{
		slv->set_ckt(*cor);
		vectord cor_res;
		slv->Update(false, update_factor, cor_res, { "GBW", "phase margin", "Av", "power" });

		//extract results and update overdesign factor

		for (Index i = 0; i < constraints.size(); i++)
		{
			if (constraints[i].name == ":gbw")
			{
				double factor = typ_res[0] / cor_res[0];
				if (factor > cor_fac[0] && isfinite<double>(factor))
				{
					cor_fac[0] = factor;
					if (factor > constraints[i].final_value)
						all_satisfied = false;
					constraints[i].overdesign_factor = factor;
					worst_res[0] = cor_res[0];
				}
			}
			else if (constraints[i].name == ":phase_m")
			{
				double factor = tan(3.14159 / 180 * (90 - cor_res[1])) / tan(3.14159 / 180 * (90 - typ_res[1]));
				if (factor > cor_fac[1] && isfinite<double>(factor))
				{
					cor_fac[1] = factor;
					if (factor > constraints[i].final_value)
						all_satisfied = false;
					constraints[i].overdesign_factor = factor;
					worst_res[1] = cor_res[1];
				}
			}
			else if (constraints[i].name == ":av")
			{
				double factor = typ_res[2] / cor_res[2];
				if (factor > cor_fac[2] && isfinite<double>(factor))
				{
					cor_fac[2] = factor;
					if (factor > constraints[i].final_value)
						all_satisfied = false;
					constraints[i].overdesign_factor = factor;
					worst_res[2] = cor_res[2];
				}
			}
			else if (constraints[i].name == ":power")
			{
				double factor = cor_res[3] / typ_res[3];
				if (factor > cor_fac[3] && isfinite<double>(factor))
				{
					cor_fac[3] = factor;
					if (factor > constraints[i].final_value)
						all_satisfied = false;
					constraints[i].overdesign_factor = factor;
					worst_res[3] = cor_res[3];
				}
			}
		}

		//If the goal is min. power. Just get the worst power value
		if (cor_res[3] > worst_res[3])
		{
			worst_res[3] = cor_res[3];
		}

	}
	slv->set_ckt(CDB);
	typ_res = worst_res;
	return all_satisfied;
}

double posterm::eval(vectord &x)
{
	double ret = coef;
	//Loop over all (super) monomial factors <=> loop over involved devices
	for (int i = 0; i < devices.size(); i++)
	{
		Device_ext* d = devices[i];
		monomial *reqMonomial = getMonomial_Model(d, i);
		//Extract indep of current device 
		vectord ev;
		for (int l = 0; l < reqMonomial->powers.size(); l++)
			ev.push_back(x[*(d->px_indices[l])] * d->x_factors[l]);
		ret *= pow(reqMonomial->leval(ev), powers[i]);
	}
	return ret;
}

vectord posterm::eval_grad(vectord &x)
{
	vectord ret(x.size());
	int i = 0;	//index of variable to derive with respect to it
	for (auto xi = x.begin(); xi != x.end(); ++xi, i++)
	{
		//get partial derivative w.r.t each variable #i
		double pd = coef;
		vector<int> found;
		vector<int> foundl;
		for (int j = 0; j < devices.size(); j++)
		{
			Device_ext* d = devices[j];
			monomial *reqMonomial = getMonomial_Model(d, j);
			//Extract indep variables location of current device 
			for (int l = 0; l < reqMonomial->powers.size(); l++)
				if (*(d->px_indices[l]) == i)
				{
					//we are in the factor in which we want to derive w.r.t one of its components
					found.push_back(j);
					foundl.push_back(l);
				}
		}
		double pdd = 0;
		int  j = 0;
		for (auto f = found.begin(); f != found.end(); f++, j++)
		{
			Device_ext* d = devices[*f];
			monomial *reqMonomial = getMonomial_Model(d, *f);
			//Extract indep location of current device 
			auto indep_index = foundl[j];
			vectord ev;
			for (int l = 0; l < reqMonomial->powers.size(); l++)
				ev.push_back(x[*(d->px_indices[l])] * d->x_factors[l]);
			//derivate the factor
			double pdtemp = powers[*f] * pow(reqMonomial->leval(ev), powers[*f] - 1) * reqMonomial->leval_grad(ev.begin(), ev.end(), indep_index);
			//Put remaining factors as they are
			int j2 = 0;
			for (auto f2 = found.begin(); f2 != found.end(); f2++, j2++)
			{
				if (f == f2)
					continue;
				Device_ext* d2 = devices[*f2];
				monomial *reqMonomial2 = getMonomial_Model(d2, *f2);

				auto indep_index2 = foundl[j2];
				vectord ev;
				for (int l = 0; l < reqMonomial->powers.size(); l++)
					ev.push_back(x[*(d2->px_indices[l])]);
				pdtemp *= pow(reqMonomial2->leval(ev), powers[*f2]);
			}
			pdd += pdtemp;
		}
		pd *= pdd;
		//Regard only factors in current posterm (all other factors are zero)
		//j is super monomial factor index
		if (found.size() > 0)
		{
			for (int j = 0; j < devices.size(); j++)
			{
				if (find(found.begin(), found.end(), j) - found.end() != 0)				//Regard the remaining factors as constant
					continue;
				Device_ext* d = devices[j];
				monomial *reqMonomial = getMonomial_Model(d, j);
				vectord ev;
				for (int l = 0; l < reqMonomial->powers.size(); l++)
					ev.push_back(x[*(d->px_indices[l])] * d->x_factors[l]);
				pd *= pow(reqMonomial->leval(ev), powers[j]);
			}
		}
		ret[i] = pd;
	}
	return ret;
}

vectord posterm::eval_hess(vectord &x)
{
	vectord ret;
	int i = 0;	//index of first variable to derive with respect to
	for (auto xi = x.begin(); xi != x.end(); ++xi, i++)
	{
		int k = 0;	//index of second variable to derive with respect to
					//Calculate only the lower triangle of Hessian matrix
		for (auto xk = x.begin(); xk != x.end() && k <= i; ++xk, k++)
		{
			//get partial derivative w.r.t each variable #i and each variable #k
			double pd = coef;
			int j;
			vector<int> foundi;		//array of "devices" (super monomial factors) indices in which we want to derive w.r.t one of its variables (var #i is one of its variables)
			vector<int> foundil;
			vector<int> foundk;		//array of "devices" (super monomial factors) indices in which we want to derive w.r.t one of its variables (var #k is one of its variables)
			vector<int> foundkl;
			//Fill foundi and foundk
			for (j = 0; j < devices.size(); j++)
			{
				Device_ext* d = devices[j];
				monomial *reqMonomial = getMonomial_Model(d, j);
				//Extract indep locations of current device for each variable i, k 
				for (int l = 0; l < reqMonomial->powers.size(); l++)
				{
					if (*(d->px_indices[l]) == i)
					{
						//we are in the factor in which we want to derive w.r.t one of its components
						foundi.push_back(j);
						foundil.push_back(l);
					}
					if (*(d->px_indices[l]) == k)
					{
						//we are in the factor in which we want to derive w.r.t one of its components
						foundk.push_back(j);
						foundkl.push_back(l);
					}
				}
			}
			vector<int> foundall = foundi;
			vector<int> foundalll = foundil;
			//foundall.insert(foundall.end(), foundk.begin(), foundk.end());	
			//concatenate foundi and foundk into foundall but with no repetitions
			int jall = -1;
			for (auto fk : foundk)
			{
				jall++;
				if (find(foundi.begin(), foundi.end(), fk) - foundi.end() != 0)
					continue;
				foundall.push_back(fk);
				foundalll.push_back(foundkl[jall]);
			}
			double pdd = 0;	//includes all terms of found factors derivatives 
							//cross sum of all possible combination of factors to be derived or not derived
			int ji = 0;
			for (auto fi = foundi.begin(); fi != foundi.end(); fi++, ji++)
			{
				int jk = 0;
				for (auto fk = foundk.begin(); fk != foundk.end(); fk++, jk++)
				{
					//performance enhancements:
					//pdtemp at fi=a and fk=b is the same as that of fi=b and fk=a
					//store evaluation of a monomial, its derivative and second derivative w.r.t some var. they are reevaluated alot
					double pdtemp;
					Device_ext* d_i = devices[*fi];
					monomial *reqMonomial_i = getMonomial_Model(d_i, *fi);

					auto indep_index_i = foundil[ji];
					vectord ev_i;
					for (int l = 0; l < reqMonomial_i->powers.size(); l++)
						ev_i.push_back(x[*(d_i->px_indices[l])] * d_i->x_factors[l]);

					if (*fi == *fk)	//both fi and fk refer to the same factor derive it twice
						if (i == k)	//derivation to the same factor and the same variable twice 
							pdtemp = powers[*fi] * (powers[*fi] - 1) * pow(reqMonomial_i->leval(ev_i), powers[*fi] - 2) * reqMonomial_i->leval_grad(ev_i.begin(), ev_i.end(), indep_index_i) * reqMonomial_i->leval_grad(ev_i.begin(), ev_i.end(), indep_index_i)
							+ powers[*fi] * pow(reqMonomial_i->leval(ev_i), powers[*fi] - 1) * reqMonomial_i->leval_2ndder(ev_i.begin(), ev_i.end(), indep_index_i);
						else
						{
							auto indep_index_k = foundkl[jk];
							pdtemp = powers[*fi] * (powers[*fi] - 1) * pow(reqMonomial_i->leval(ev_i), powers[*fi] - 2) * reqMonomial_i->leval_grad(ev_i.begin(), ev_i.end(), indep_index_i) * reqMonomial_i->leval_grad(ev_i.begin(), ev_i.end(), indep_index_k)
								+ powers[*fi] * pow(reqMonomial_i->leval(ev_i), powers[*fi] - 1) * reqMonomial_i->leval_2ndder(ev_i.begin(), ev_i.end(), indep_index_i, indep_index_k);
						}
					else
					{
						Device_ext* d_k = devices[*fk];
						monomial *reqMonomial_k = getMonomial_Model(d_k, *fk);

						auto indep_index_k = foundkl[jk];
						vectord ev_k;
						for (int l = 0; l < reqMonomial_k->powers.size(); l++)
							ev_k.push_back(x[*(d_k->px_indices[l])] * d_k->x_factors[l]);

						pdtemp = powers[*fi] * pow(reqMonomial_i->leval(ev_i), powers[*fi] - 1) * reqMonomial_i->leval_grad(ev_i.begin(), ev_i.end(), indep_index_i)
							* powers[*fk] * pow(reqMonomial_k->leval(ev_k), powers[*fk] - 1) * reqMonomial_k->leval_grad(ev_k.begin(), ev_k.end(), indep_index_k);
					}
					//put remaining factors in found but not fi or fk which we have already derived and have been put
					int j2 = 0;
					for (auto f2 = foundall.begin(); f2 != foundall.end(); f2++, j2++)
					{
						if (*fi == *f2 || *fk == *f2)
							continue;
						Device_ext* d2 = devices[*f2];
						monomial *reqMonomial2 = getMonomial_Model(d2, *f2);

						auto indep_index2 = foundalll[j2];
						vectord ev2;
						for (int l = 0; l < reqMonomial2->powers.size(); l++)
							ev2.push_back(x[*(d2->px_indices[l])] * d2->x_factors[l]);
						pdtemp *= pow(reqMonomial2->leval(ev2), powers[*f2]);
					}
					pdd += pdtemp;
				}
			}
			pd *= pdd;
			if (foundi.size() > 0 || foundk.size() > 0)
			{
				//put remaining factors not in foundi nor in foundk
				for (int l = 0; l < devices.size(); l++)
				{
					if (find(foundi.begin(), foundi.end(), l) - foundi.end() != 0 || find(foundk.begin(), foundk.end(), l) - foundk.end() != 0)
						continue;
					Device_ext* d = devices[l];
					monomial *reqMonomial = getMonomial_Model(d, l);
					//Extract indep location of current device 
					vectord ev;
					for (int l = 0; l < reqMonomial->powers.size(); l++)
						ev.push_back(x[*(d->px_indices[l])] * d->x_factors[l]);
					//Regard the remaining factors as constant
					pd *= pow(reqMonomial->leval(ev), powers[l]);
				}
			}
			ret.push_back(pd);
		}
	}
	return ret;
}

double posynomial::eval(vectord &x)
{
#ifdef MEAS_TIME
	auto start = chrono::high_resolution_clock::now();
#endif
	double sum = 0;
	for (auto m : terms)
		sum += m.eval(x);
	sum *= overdesign_factor;
	final_value = sum;
#ifdef MEAS_TIME
	auto stop = chrono::high_resolution_clock::now();
	eval_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif
	return sum;
}

vectord posynomial::eval_grad(vectord &x)
{
#ifdef MEAS_TIME
	auto start = chrono::high_resolution_clock::now();
#endif
	vectord ret(x.size());
	for (auto m : terms)
	{
		vectord tmp = m.eval_grad(x);
		auto ti = tmp.begin();
		for (auto ri = ret.begin(); ri != ret.end(); ri++, ti++)
			*ri += *ti;
	}
	ret = ret * overdesign_factor;
#ifdef MEAS_TIME
	auto stop = chrono::high_resolution_clock::now();
	grad_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif
	return ret;
}

vectord posynomial::eval_hess(vectord &x, double multiplier = 1)
{
#ifdef MEAS_TIME
	auto start = chrono::high_resolution_clock::now();
#endif
	//for w.r.t 6, 7 
	//value is 3.9827e+8
	//value should be 1.561e+7
	//x = {1e-5, 1e-6, 2e-6, 1e-6, 2e-6, 1e-6, 2e-6, 1e-6};
	size_t n = x.size();
	vectord ret(n * (n + 1) / 2);
	for (auto m : terms)
	{
		vectord tmp = m.eval_hess(x);
		auto ti = tmp.begin();
		for (auto ri = ret.begin(); ri != ret.end(); ri++, ti++)
			*ri += *ti;
	}
	if (multiplier != 1)
		for (auto ri = ret.begin(); ri != ret.end(); ri++)
			*ri *= multiplier;
	ret = ret * overdesign_factor;
#ifdef MEAS_TIME
	auto stop = chrono::high_resolution_clock::now();
	hess_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif
	return ret;
}

vectord posynomial::eval_hess(vectord &x, const vectord &add, double multiplier)
{
#ifdef MEAS_TIME
	auto start = chrono::high_resolution_clock::now();
#endif
	size_t n = x.size();
	vectord ret(n * (n + 1) / 2);
	for (auto m : terms)
	{
		vectord tmp = m.eval_hess(x);
		auto ti = tmp.begin();
		for (auto ri = ret.begin(); ri != ret.end(); ri++, ti++)
			*ri += *ti;
	}
	auto addi = add.begin();
	for (auto ri = ret.begin(); ri != ret.end(); ri++, addi++)
		*ri = *ri * multiplier * overdesign_factor + *addi;
#ifdef MEAS_TIME
	auto stop = chrono::high_resolution_clock::now();
	hess_clock += chrono::duration_cast<chrono::microseconds>(stop - start);
#endif
	return ret;
}

void GPOptimizer::run_test(int maxIter)
{


}
