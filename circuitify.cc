#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gmp.h>
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include <stdexcept>
#include <gmpxx.h>
#include <iostream>
using namespace std;

#define MODULUS 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
#define HALFMODULUS = ((MODULUS + 1)/2)
#define COST_SCALAR_MUL 5
#define COST_SCALAR_NEG 2
#define COST_SCALAR_COPY 1

typedef tuple<char, int> var_tuple;

vector<struct mul> mul_data;

int mul_count = 0;
int temp_count = 0;


class Vars {
public:
	map<int, mpz_class> var_map;

	void add_var(char type, int idx, int val) {
		int s;
		if (type == 'L')
			s = 0;
		else if (type == 'R') 
			s = 1;
		else if (type == 'O')
			s = 2;
		else if (type == 'T')
			s = 3;
		else
			throw invalid_argument("Invalid variable type");
		var_map[4*idx + s] = val;
	}

	void add(Vars& other) {
		for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				var_map[x.first] = x.second;
			} else {
				var_map[x.first] += x.second;
			}
		}
	}

	void sub(Vars& other) {
		for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				var_map[x.first] = x.second;
			} else {
				var_map[x.first] -= x.second;
			}
		}
	}

	void mul(mpz_class v) {
		for (auto const&x : var_map) {
			var_map[x.first] *= v;
		}
	}

	void div(mpz_class v) {
		for (auto const&x : var_map) {
			var_map[x.first] /= v;
		}
	}
};

class Linear {

	public:
	mpz_class val;
	bool is_const;
	Vars vars; 

	/*Linear(int v) {
		val = v;
	}*/

	void add_var(char type, int idx, int val) {
		vars.add_var(type, idx, val);
	}

	void add(Linear& other) {
		val += other.val;
		vars.add(other.vars);
	}

	void sub(Linear& other) {
		val -= other.val;
		vars.sub(other.vars);
	}

	void mul(mpz_class v) {
		val *= v;
		vars.mul(v);
	}

	void div(mpz_class v) {
		val /= v;
		vars.div(v);
	}
};

vector<Linear> eqs;

struct mul {
	mpz_class l;
	mpz_class r;
	mpz_class o;
};

struct multiplication {
	Linear l;
	Linear r;
	Linear o;
};

mpz_class modinv(mpz_class n) {
	
	return 5;//(int)pow(n, MODULUS-2);
}

void new_mul(mpz_class l, mpz_class r, Linear& nl, Linear& nr, Linear& no)
{
	mpz_class o = l*r;
	struct mul m = {l, r, o};
	mul_data.push_back(m);
	nl.val = l;
	nl.is_const = false;
	nl.add_var('L', mul_count, 1);
	nr.val = r;
	nr.is_const = false;
	nr.add_var('R', mul_count, 1);
	no.val = o;
	no.is_const = false;
	no.add_var('O', mul_count, 1);
	mul_count += 1;
}

void new_temp(mpz_class v, Linear& nt) {
	nt.val = v;
	nt.is_const = false;
	nt.add_var('T', mul_count, 1);
	temp_count += 1;
}

void new_const(mpz_class v, Linear& nc) {
	nc.val = v;
	nc.is_const = true;
}

// currently mutates l and/or r
Linear new_multiplication(Linear& l, Linear& r, bool addeqs = true) {
	if (l.is_const) {
		r.mul(l.val);
		return r;
	}
	if (r.is_const) {
		l.mul(r.val);
		return l;
	}
	if (r.val < l.val) {
		Linear tmp = l;
		l = r;
		r = tmp;
	}
	Linear lv = Linear();
	Linear rv = Linear();
	Linear ret = Linear();
	new_mul(l.val, r.val, lv, rv, ret);
	l.sub(lv);
	eqs.push_back(l);
	if (addeqs){
		r.sub(rv);
		eqs.push_back(r);
	}
	return ret;
}

Linear new_division(Linear& l, Linear& r) {
	if (r.is_const) {
		l.div(r.val);
		return l;
	}
	Linear lv = Linear();
	Linear rv = Linear();
	Linear ret = Linear();
	new_mul(l.val * modinv(r.val), r.val, lv, rv, ret);
	l.sub(lv);
	r.sub(rv);
	eqs.push_back(l);
	eqs.push_back(r);
	return ret;
}


int main() {
//  printf("%d\n", modinv(5));
  
  mpz_class n;
  n = 0;
  printf("n = ");
  cout << n;
  printf("\n");
  return 0;
}
