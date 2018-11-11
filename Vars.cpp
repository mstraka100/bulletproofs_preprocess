#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <sstream>
#include <gmpxx.h>
#include "utils.h"
#include "Vars.h"

using namespace std;


#define COST_SCALAR_MUL 5
#define COST_SCALAR_NEG 2
#define COST_SCALAR_COPY 1

	// deletes variables with value 0
	void Vars::reduce() {
		for(map<int, mpz_class>::iterator it = var_map.begin(); it != var_map.end();) {
			if((it->second % mod) == 0) {
				it = var_map.erase(it);
			} else {
				it++;
			}
		}
	}

	void Vars::add_var(char type, int idx, int val) {
		int s = var_offset(type);
		var_map[4*idx + s] = val % mod;
	}

	bool Vars::has_var(char type, int idx) {
		int s = var_offset(type);
		int key = 4*idx + s;
		return var_map.count(key) != 0;
	}

	mpz_class Vars::get_var(char type, int idx) {
		int s = var_offset(type);
		int key = 4*idx + s;
		return var_map[key];
	}

	void Vars::index_temp_vars(map<int, vector<int>>& index, int pos) {
		for (auto const&x : var_map) {
			if (var_type(x.first) == 'T') {
				int idx = var_idx(x.first);
				if (index.count(idx) == 0)
					index[idx] = {pos};
				else 
					index[idx].push_back(pos);
			}
		}
	}

	int Vars::num_vars() {
		return var_map.size();
	}

	void Vars::add(Vars& other) {
		for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				var_map[x.first] = x.second % mod;
			} else {
				var_map[x.first] = (var_map[x.first] + x.second) % mod;
			}
		}
		reduce();
	}

	void Vars::sub(Vars& other) {
		for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				var_map[x.first] = (mod - x.second) % mod;
			} else {
				var_map[x.first] = (var_map[x.first] - x.second + mod) % mod;
			}
		}
		reduce();
	}

	void Vars::mul(mpz_class v) {
		for (auto const&x : var_map) {
			var_map[x.first] = (var_map[x.first] * v) % mod;
		}
	}

	void Vars::div(mpz_class v) {
		for (auto const&x : var_map) {
			var_map[x.first] = (var_map[x.first] / v) % mod;
		}
	}	

	bool Vars::is_empty() {
		return var_map.empty();
	}

	bool Vars::is_zero() {
		for (auto const&x : var_map) {
			if (x.second != 0) {
				return false;
			}
		}
		return true;
	}

	void Vars::str(mpz_class constant) {
		vector<string> terms;
		for (auto const&e : var_map) {
			string name = var_str(e.first);
			mpz_class val = e.second;
			stringstream s;
			if (val == 1)
				terms.push_back(name);
			else if ((val + 1) % mod == 0)
				terms.push_back("-" + name);
			else if (mod - val < 10000) {
				s << "-" << mod-val << "*" << name;
				terms.push_back(s.str());
			} else {
				s << (val) << "*" << name;
				terms.push_back(s.str());
			}
		}
		stringstream s;
	    if (terms.empty() || constant != 0) {
			if (mod - constant < 10000) {
				s << "-" << mod-constant;
				terms.push_back(s.str());
			}
			else {
				s << constant;
				terms.push_back(s.str());
			}
		}
		
		for (auto const&x : terms) {
			cout << x << " + ";
		}
	}

	int Vars::cost() {
		int total = 0;
		map<mpz_class, int> negs;
		map<mpz_class, int> counts;
		mpz_class high = 0;
		for (auto& e : var_map) {
			mpz_class val = e.second;
			total += 1;
			int negate = 0;
			mpz_class diff = mod - val;
			if (diff < val) {
				val = diff;
				negate = 1;
			}
			if (counts.count(val) != 0) {
				counts[val] += 1;
				negs[val] += negate;
			} else {
				counts[val] = 1;
				negs[val] = negate;
			}
			if (high == 0 || counts[val] > counts[high]) {
				high = val;
			}
		}
		int plus_one = max(negs[high], counts[high] - negs[high]);
		int neg_one = counts[high] - plus_one;
		return COST_SCALAR_MUL * (total - plus_one - neg_one) + COST_SCALAR_NEG * neg_one + COST_SCALAR_COPY;
	}

	map<int, mpz_class>::iterator Vars::vars_begin() {
		return var_map.begin();
	}

	map<int, mpz_class>::iterator Vars::vars_end() {
		return var_map.end();
	}