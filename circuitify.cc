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
#include <array>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
using namespace std;

#define COST_SCALAR_MUL 5
#define COST_SCALAR_NEG 2
#define COST_SCALAR_COPY 1

typedef tuple<char, int> var_tuple;

vector<struct mul> mul_data;


int mul_count = 0;
int temp_count = 0;
int bit_count = 0;
mpz_class mod = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141_mpz;

regex var_re = regex("[A-Za-z_][0-9a-zA-Z_]*");
regex secret_re = regex("#(-?[0-9]+)");
regex num_re = regex("[0-9]+");

class Vars {

private:
	int var_offset(char type) {
		if (type == 'L')
			return 0;
		else if (type == 'R') 
			return 1;
		else if (type == 'O')
			return 2;
		else if (type == 'T')
			return 3;
		else
			throw invalid_argument("Invalid variable type");
	}

public:
	map<int, mpz_class> var_map;


	void add_var(char type, int idx, int val) {
		int s = var_offset(type);
		var_map[4*idx + s] = val;
	}

	bool has_var(char type, int idx) {
		int s = var_offset(type);
		int key = 4*idx + s;
		return var_map.count(key) != 0;
	}

	mpz_class get_var(char type, int idx) {
		int s = var_offset(type);
		int key = 4*idx + s;
		return var_map[key];
	}

	int num_vars() {
		return var_map.size();
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

	void to_str() {
		for (std::map<int,mpz_class>::iterator it=var_map.begin(); it!=var_map.end(); ++it)
    		std::cout << it->first << " => " << it->second << '\n';
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

	bool has_var(char type, int idx) {
		return vars.has_var(type, idx);
	}

	int num_vars() {
		return vars.num_vars();
	}

	mpz_class get_var(char type, int idx) {
		return vars.get_var(type, idx);
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

	void to_str() {
		vars.to_str();
	}
};

vector<Linear> eqs;
map<string, Linear> varset;

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

mpz_class modinv(mpz_class x) {
	mpz_t ret;
	mpz_init(ret);
	mpz_class n = mod - 2;
	mpz_powm(ret, x.get_mpz_t(), n.get_mpz_t(), mod.get_mpz_t());
	return mpz_class(ret);
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

// mutates l and/or r
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

Linear new_xor(Linear& l, Linear& r) {
	Linear lv = Linear();
	Linear rv = Linear();
	Linear mul = Linear();
	new_multiplication(l, r);
	l.add(r);
	mul.mul(2);
	l.sub(mul);
	return l;
}

string clean_expr(string s) {
	// TODO: strip s
	boost::algorithm::trim(s);
	if (s == "" || s[0] != '(' || s.back() == ')')
		return s;
	int depth = 1;
	for (int i = 1; i < s.length()-1; i++) {
		if (s[i] == '(') {
			depth += 1;
		} else if (s[i] == ')') {
			depth -= 1;
			if (depth == 0)
				return s;
		}
	}
	return clean_expr(s.substr(1, s.length()-2));
}

struct expr {
	string l;
	string op;
	string r;
};

bool split_expr_binary(string& s, vector<string> ops, struct expr& result) {
	//int i = 1;
	int depth = 0;
	for (int i = s.length(); i > 1; i--) {
		if (s[i-1] == ')')
			depth += 1;
		else if (s[i-1] == '(')
			depth -= 1;
		else if (depth == 0){
			for (const string& op: ops) {
				if (i-op.length() >= 1 && s.substr(i-op.length(), op.length()) == op) {
					result.l = clean_expr(s.substr(0, i-op.length()));
					result.r = clean_expr(s.substr(i));
					result.op = op;
					return true;			
				}
			}
		}
	}
	return false;
}

Linear parse_expression(string s) {
	cout << s << endl;
	s = clean_expr(s);
	bool split;
	if (s == "")
		throw invalid_argument("Empty expresion");
	struct expr sp;
	vector<string> delim = {"^"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l);
		Linear right = parse_expression(sp.r);
		Linear ret = new_xor(left, right);
		// TODO: assert result correct
		return ret;
	}
	delim = {"+", "-"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l);
		Linear right = parse_expression(sp.r);
		if (sp.op == "+")
			left.add(right);
		else
			left.sub(right);
		return left;
	}
	delim = {"*", "/"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l);
		Linear right = parse_expression(sp.r);
		if (sp.op == "*") {
			Linear ret = new_multiplication(left, right);
			return ret;
		} else {
			Linear ret = new_division(left, right);
			return ret;
		}
	}
	if (s.length() > 5 && s.substr(0,5) == "bool(") {
		Linear ret = parse_expression(s.substr(4));
		Linear tmp = ret;
		Linear tmp2 = ret;
		Linear one;
		new_const(1, one);
		tmp2.sub(one);
		new_multiplication(tmp, tmp2, false);
		bit_count += 1;
		return ret;
	}
	if (s[0] == '-') {
		Linear ret = parse_expression(s.substr(1));
		ret.mul(mod-1);
	}
	if (regex_match(s, var_re)) {
		if (varset.count(s) != 0)
			return varset[s];
		else
			throw invalid_argument("Variable not defined");
	}
	if (regex_match(s, secret_re)) {
		Linear ret;
		int val = stoi(s.substr(1));
		mpz_class mpz_val = val;
		new_temp(mpz_val, ret);
		return ret;
	}
	if (regex_match(s, num_re)) {
		Linear ret;
		int val = stoi(s);
		mpz_class mpz_val = val;
		new_const(mpz_val, ret);
		return ret;
	}
	throw invalid_argument("Cannot parse expression");
}

vector<Linear> parse_expressions(string s) {
	struct expr sp;
	vector<string> delim = {","};
	bool split = split_expr_binary(s, delim, sp);
	if (split) {
		vector<Linear> l = parse_expressions(sp.l);
		Linear r = parse_expression(sp.r);
		l.push_back(r);
		return l;
	}
	Linear ret_val = parse_expression(s);
	vector<Linear> ret;
	ret.push_back(ret_val);
	return ret;
}


//void split_by_delimiter(constconst string delim, )

bool all_const(vector<Linear> v) {
	for (int i = 0; i < v.size(); i++) {
		if (!v[i].is_const)
			return false;
	}
	return true;
}

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
void parse_statement(string& s) {
	//TODO: strip s
	if (s.length() > 6 && s.substr(0,6) == "debug "){
		throw invalid_argument("TODO: implement debug");
	}
	vector<string> assn_ops = {":=", "=:", "==", "="};
	struct expr sp;
	bool match = split_expr_binary(s, assn_ops, sp);
	if (match) {
		string left = sp.l;
		string op = sp.op;
		string right = sp.r;

		if (op == ":=") {
			vector<string> bits;
			//TODO: Abstract out splitting string
			boost::char_separator<char> sep{","};
  			tokenizer tok{left, sep};
			for (const auto &b : tok)
				bits.push_back(b);
			//TODO: Check bit names are valid
			Linear val = parse_expression(right);
			//TODO: assert bits length <= 256
			vector<Linear> bitvars;
			if (val.is_const) {
				for (int i = 0;i < bits.size(); i++) {
					Linear l;
					mpz_class v = (val.val >> i) & 1;
					new_const(v,l);
					bitvars.push_back(l);
				} 
			} else {
				for (int i = 0; i < bits.size(); i++) {
					Linear l;
					mpz_class v = (val.val >> i) & 1;
					new_temp(v,l);
					bitvars.push_back(l);
					Linear tmp = l;
					Linear tmp2 = l;
					Linear one;
					new_const(1, one);
					tmp2.sub(one);
					Linear eq = new_multiplication(tmp, tmp2);
					eqs.push_back(eq);
					tmp = l;
					tmp.mul(1 << i);
					val.sub(tmp);
				}
			}
			for (int i = 0; i < bits.size(); i++) {
				varset[bits[i]] = bitvars[i];
			}
		} else if (op == "=:") {
			vector<Linear> bits;
			bits = parse_expressions(right);
			// TODO: Check varible name on left valid, len of bits <=256
			mpz_class real = 0;
			for (int i = 0; i < bits.size(); i++) {
				//TODO: Assert bit either 0 or 1
				real = real + bits[i].val * (1 << i);
				if (all_const(bits)) {
					Linear val;
					new_const(real, val);
					varset[left] = val;
				} else {
					Linear val;
					new_temp(real, val);
					varset[left] = val;
					for (int i = 0; i < bits.size(); i++) {
						Linear tmp = bits[i];
						tmp.mul(1 << i);
						val.sub(tmp);
					}
					eqs.push_back(val);
				}
			}
		} else if (op == "=") {
			//TODO: check left is proper variable name
			Linear ex = parse_expression(right);
			varset[left] = ex;
		} else if (op == "==") {
			Linear l = parse_expression(left);
			Linear r = parse_expression(right);
			l.sub(r);
			//TODO: check l.val == r.val
			eqs.push_back(l);
		}
	}
}


void pivot_variable(char type, int idx, bool eliminate = false) {
	int c = 0;
	int cc = 0;
	int low = -1;
	Linear leq;
	vector<int> vec;
	for (int i = 0; i < eqs.size(); i++) {
		if (eqs[i].has_var(type, idx)) {
			vec.push_back(i);
			cc += 1;
			if (low == -1 || c > eqs[i].num_vars()) {
				low = i;
				c = eqs[i].num_vars();
				Linear tmp = eqs[i];
				mpz_class v = eqs[i].get_var(type, idx);
				mpz_class inv = modinv(v);
				tmp.mul(inv);
				leq = eqs[i];
				break;
			}
		}
	}
	if (cc > 1) {
		for (auto i: vec) {
			if (i != low) {
				Linear tmp = leq;
				tmp.mul(eqs[i].get_var(type, idx));
				eqs[i].sub(tmp);
			}
		}
	}
	if (eliminate and cc > 0) {
		eqs.erase(eqs.begin()+low);
	}
}

void eliminate_temps() {
	for (int i = 0; i < temp_count; i++) {
		pivot_variable('T', i, true);
	}
}


int main() {
//  printf("%d\n", modinv(5));
  
/*  mpz_class n;
  n = 0;
  printf("n = ");
  cout << n;
  printf("\n");*/

	vector<string> input = {"v1 = #2", "v2 = #5", "v3 = v1 * v2", "v4 = v3 * v3"};

	for (int i = 0; i < input.size(); i++)
		parse_statement(input[i]);

	//printf("%d multiplications, %d temporaries, %lu constraints", mul_count, temp_count, eqs.size());
	
	printf("%d multiplications, %d temporaries, %lu constraints", mul_count, temp_count, eqs.size());
/*	for (int i = 0; i < eqs.size(); i++) {
		eqs[i].to_str();
	} */

	eliminate_temps();

	printf("did it!");

 	return 0;
}
