#include "Vars.h"
#include "Linear.h"
#include "utils.h"


	void Linear::add_var(char type, int idx, int val) {
		vars.add_var(type, idx, val);
	}

	bool Linear::has_var(char type, int idx) {
		return vars.has_var(type, idx);
	}

	int Linear::num_vars() {
		return vars.num_vars();
	}

	mpz_class Linear::get_var(char type, int idx) {
		return vars.get_var(type, idx);
	}

	void Linear::add(Linear& other) {
		real = (real + other.real) % mod;
		constant = (constant + other.constant) % mod;
		vars.add(other.vars);
	}

	void Linear::sub(Linear& other) {
		real = (real - other.real + mod) % mod;
		constant = (constant - other.constant + mod) % mod;
		vars.sub(other.vars);
	}

	void Linear::mul(mpz_class v) {
		real = (real * v) % mod;
		constant = (constant * v) % mod;
		vars.mul(v);
	}

	void Linear::div(mpz_class v) {
		v = modinv(v, mod);
		mul(v);
	}

	void Linear::index_temp_vars(map<int, vector<Linear*>>& index) {
		vars.index_temp_vars(index, *this);
	}

	void Linear::assign_temp_vars(Linear& other, map<int, vector<Linear*>>& index) {
		other.vars.index_temp_vars(index, *this, true);
	}

	void Linear::to_str() {
		vars.str(constant);
	}

	int Linear::equation_cost() {
		return vars.cost();
	}

	bool Linear::is_const() {
		return vars.is_empty();
	}

	bool Linear::is_zero() {
		if (constant != 0) {
			return false;
		}
		return vars.is_zero();
	}

	bool Linear::has_temp_var() {
		return vars.has_temp_var();
	}

	bool Linear::has_several_temps() {
		return vars.has_several_temps();
	}

	map<int, mpz_class>::iterator Linear::vars_begin() {
		return vars.vars_begin();
	}

	map<int, mpz_class>::iterator Linear::vars_end() {
		return vars.vars_end();
	}