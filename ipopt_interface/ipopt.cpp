#include <iostream>
#include "model.h"
#include "modeldb.h"
#include "circuit.h"
#include "optimizergp.h"
#include "ipopt.h"

#define logarithmic

//IPMOpt

bool IPMOpt::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	n = (Index)(this->vars->size());
	m = (Index)(this->constrs->size());

	// in this example the jacobian is dense and contains {} nonzeros
	//nnz_jac_g = n * m;
	vector<Number> xvec(n);
	for (auto xi = xvec.begin(); xi != xvec.end(); xi++)
		*xi = 10;
	vector<Number> tres;
	int nnz = 0;
	for (auto p : *constrs)
	{
		tres = p.eval_grad(xvec);
		//append to end result (total size is n*m)
		for (auto ri = tres.begin(); ri != tres.end(); ri++)
		{
			if (*ri != 0)
				nnz++;
		}
	}
	nnz_jac_g = nnz;
	// the hessian is also dense and has {n*n} total nonzeros, but we
	// only need the lower left corner (since it is symmetric)
	nnz_h_lag = n * (n + 1) / 2;

	// use the C style indexing (0-based)
	index_style = TNLP::C_STYLE;

	return true;
}

// returns the variable bounds
bool IPMOpt::get_bounds_info(Index n, Number* x_l, Number* x_u,
	Index m, Number* g_l, Number* g_u)
{
	// here, the n and m we gave IPOPT in get_nlp_info are passed back to us.

	// the variables have lower and upper boundss
	for (Index i = 0; i < n; i++) {
		x_l[i] = (*vrngs)[i].begin;
		x_u[i] = (*vrngs)[i].end;
	}

	// set constraints bounds
	for (Index i = 0; i < m; i++) {
		switch ((*constrs_types)[i])
		{
			case -1:
				g_l[i] = 0;
				g_u[i] = 1;
				break;
			case 0:
				g_l[i] = 1;
				g_u[i] = 1;
				break;
			case 1:
				g_l[i] = 1;
				g_u[i] = 0;
				break;
		default:
			break;
		}
	}
	return true;
}

// returns the initial point for the problem
bool IPMOpt::get_starting_point(Index n, bool init_x, Number* x,
	bool init_z, Number* z_L, Number* z_U,
	Index m, bool init_lambda,
	Number* lambda)
{
	// Here, we assume we only have starting values for x
	//assert(init_x == true);
	//assert(init_z == false);
	//assert(init_lambda == false);

	// initialize to the given starting point
	for (Index i = 0; i < n; i++) {
		x[i] = (*vars)[i];
	}

	return true;
}


// returns the value of the objective function
bool IPMOpt::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	raw_vars = x;
	vector<Number> xvec(x, x + n);
	obj_value = goal->eval(xvec);
	return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IPMOpt::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	vector<Number> xvec(x, x + n);
	vector<Number> res = goal->eval_grad(xvec);
	int i = 0;
	for (auto ri = res.begin(); ri != res.end(); ri++, i++)
		grad_f[i] = *ri;
	return true;
}

// return the value of the constraints: g(x)
bool IPMOpt::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	vector<Number> xvec(x, x + n);
	int i = 0;
	for (auto p = constrs->begin(); p != constrs->end(); p++)
	{
		g[i] = (*p).eval(xvec);
		i++;
	}
	return true;
}

// return the structure or values of the jacobian
bool IPMOpt::eval_jac_g(Index n, const Number* x, bool new_x,
	Index m, Index nele_jac, Index* iRow, Index *jCol,
	Number* values)
{
	if (values == NULL) {
		// return the structure of the jacobian
		// this particular jacobian is dense
		//for (int i = 0; i < m; i++) {
		//	for (int j = 0; j < n; j++)
		//	{
		//		iRow[n*i + j] = i;
		//		jCol[n*i + j] = j;
		//	}
		//}

		vector<Number> xvec(n);
		for (auto xi = xvec.begin(); xi != xvec.end(); xi++)
			*xi = 10;
		vector<Number> tres;
		int i = 0;
		int nnz = 0;
		for (auto p = constrs->begin(); p!= constrs->end(); p++, i++)
		{
			tres = (*p).eval_grad(xvec);
			int j = 0;
			for (auto ri = tres.begin(); ri != tres.end(); ri++, j++)
			{
				if (*ri != 0)
				{
					iRow[nnz] = i;
					jCol[nnz] = j;
					nnz++;
				}
			}
		}
	}
	else {
		// return the values of the jacobian of the constraints
		vector<Number> xvec(x, x + n);
		vector<Number> tres;
		int i = 0;
		for (auto p = constrs->begin(); p != constrs->end(); p++)
		{
			tres = (*p).eval_grad(xvec);
			//append to end result (total size is n*m)
			for (auto ri = tres.begin(); ri != tres.end(); ri++)
			{
				if (*ri != 0)
				{
					values[i] = *ri;
					i++;
				}
			}
		}
	}
	return true;
}

//return the structure or values of the hessian
bool IPMOpt::eval_h(Index n, const Number* x, bool new_x,
	Number obj_factor, Index m, const Number* lambda,
	bool new_lambda, Index nele_hess, Index* iRow,
	Index* jCol, Number* values)
{
	if (values == NULL) {
		// return the structure. This is a symmetric matrix, fill the lower left
		// triangle only.

		// the hessian for this problem is actually dense
		Index idx = 0;
		for (Index row = 0; row < n; row++) {
			for (Index col = 0; col <= row; col++) {
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
			}
		}

		//assert(idx == nele_hess);
	}
	else {
		// return the values. This is a symmetric matrix, fill the lower left
		// triangle only

		// fill the objective portion
		vector<Number> xvec(x, x + n);
		auto hessian = goal->eval_hess(xvec, obj_factor);
		int j = 0;
		for (auto p : *constrs)
		{
			hessian = p.eval_hess(xvec, hessian, lambda[j]);
			j++;
		}
		int i = 0;
		for (auto ri = hessian.begin(); ri != hessian.end(); ri++, i++)
			values[i] = *ri;
	}

	return true;
}

bool IPMOpt::intermediate_callback(AlgorithmMode mode, Index iter, Number obj_value,
	Number inf_pr, Number inf_du, Number mu, Number d_norm, Number regularization_size,
	Number alpha_du, Number alpha_pr, Index ls_trials, const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
{
	//return true;
//	Ipopt::TNLPAdapter* tnlp_adapter = NULL;
	if (ip_cq == NULL)
		return true;		//restoration phase
//	else
//	{

//		Ipopt::OrigIpoptNLP* orignlp;
//		orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
//		if (orignlp != NULL)
//			tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
//	}
//	tnlp_adapter->ResortX(*ip_data->curr()->x(), primals);
//	vectord curr_vars(primals, vars->size());

//determine maximum relative change
	vectord curr_vars(vars->size());
	for (int j = 0; j < vars->size(); j++)
		curr_vars[j] = raw_vars[j];

#ifdef logarithmic
	double target = log(1 + MAX_CHANGE_PERCENT / 100 * shrink_factor);
	vectord diff = curr_vars - init_point;
#else
	double target = MAX_CHANGE_PERCENT / 100 * shrink_factor;
	vectord diff = (curr_vars - init_point) / init_point;
#endif
	for (auto e = diff.begin(); e != diff.end(); e++)
	{
		if (abs(*e) > target)
			return false;
	}
	return true;
}

void IPMOpt::finalize_solution(SolverReturn status,
	Index n, const Number* x, const Number* z_L, const Number* z_U,
	Index m, const Number* g, const Number* lambda,
	Number obj_value,
	const IpoptData* ip_data,
	IpoptCalculatedQuantities* ip_cq)
{

	// here is where we would store the solution to variables, or write to a file, etc
	// so we could use the solution.

	// For this example, we write the solution to the referenced vars
	/*
	std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
	for (Index i = 0; i < n; i++) {
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
	}
	*/
	for (Index i = 0; i < n; i++)
		(*vars)[i] = x[i];

#ifdef logarithmic
	double target = log(1 + MAX_CHANGE_PERCENT / 100 * shrink_factor);
	vectord diff = *vars - init_point;
#else
	double target = MAX_CHANGE_PERCENT / 100 * shrink_factor;
	vectord diff = (curr_vars - init_point) / init_point;
#endif

	int i = 0;
	for (auto e = diff.begin(); e != diff.end(); e++, i++)
	{
		if (abs(*e) > target)
		{
			//revert the variable to its change limit
#ifdef logarithmic
			(*vars)[i] = init_point[i] + abs(diff[i]) / diff[i] * target;
#else
			(*vars)[i] = init_point[i] * (1 + abs(diff[i]) / diff[i] * target);
#endif
		}
	}

	//reset initial point to current point
	init_point = *vars;

	/*
	std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
	for (Index i = 0; i < n; i++) {
		std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
	}
	for (Index i = 0; i < n; i++) {
		std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
	}
	*/
	std::cout << std::endl << std::endl << "Objective value:  ";
	std::cout << "f(x*) = " << obj_value << std::endl;

	std::cout << std::endl << "Final value of the constraints:" << std::endl;
	for (Index i = 0; i < m; i++) {
		std::cout << "g(" << i << ") = " << g[i] << std::endl;
	}
	std::cout << std::endl;
}