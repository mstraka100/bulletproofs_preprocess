#ifndef LINEAR_H
#define LINEAR_H

#include <map>
#include <vector>
#include "Vars.h"

using namespace std;

class Linear {

public:
	mpz_class real;
	mpz_class constant;
	Vars vars;
	bool eliminated = false;

	void add_var(char type, int idx, int val);

	bool has_var(char type, int idx);

	int num_vars();

	mpz_class get_var(char type, int idx);

	void add(Linear& other);

	void sub(Linear& other);

	void mul(mpz_class v);

	void div(mpz_class v);

	void index_temp_vars(map<int, vector<int>>& index, int idx);

	void assign_temp_vars(Linear& other, map<int, vector<int>>& index, int idx);

	void to_str();

	int equation_cost();

	bool is_const();

	bool is_zero();

	bool has_temp_var();

	bool has_several_temps();

	map<int, mpz_class>::iterator vars_begin();

	map<int, mpz_class>::iterator vars_end();
};

#endif