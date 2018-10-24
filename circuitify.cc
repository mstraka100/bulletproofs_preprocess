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
map<string, Linear> varset;

int mul_count = 0;
int temp_count = 0;
int bit_count = 0;
mpz_class mod = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141_mpz;

regex var_re = regex('[A-Za-z_][0-9a-zA-Z_]*');
regex secret_re = regex('#(-?[0-9]+)');
regex num_re = regex('[0-9]+');

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

mpz_class modinv(mpz_class x) {
	mpz_t ret;
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
	s = clean_expr(s);
	bool split;
	if (s == "")
		throw invalid_argument("Empty expresion");
	struct expr sp;
	vector<string> delim = {"^"}
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l);
		Linear right = parse_expression(sp.r);
		Linear ret = new_xor(left, right)
		// TODO: assert result correct
		return ret;
	}
	delim = {"+", "-"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l);
		Linear right = parse_expression(sp.r);
		if (sp.op == "+")
			sp.l.add(sp.r);
		else
			sp.l.sub(sp.r)
		return sp.l;
	}
	delim = {"*", "/"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l);
		Linear right = parse_expression(sp.r);
		if (op == "*") {
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
		temp2.sub(one);
		new_multiplication(tmp, tmp2, false);
		bit_count += 1;
		return ret;
	}
	if (s[0] == '-') {
		Linear ret = parse_expression(s.substr(1))
		ret.mul(mod-1);
	}
	if (regex_match(s, var_re)) {
		if (varset.count(s) != 0)
			return varset[s];
		else
			throw invalid_argument("Variable not defined");
	}
	
}

/*void parse_expressions(string s) {
	sruct expr sp;
	vector<string> delim = {","};
	bool split = split_expr_binary(s, delim);
	if (split) {
		return parse_expressions
	}
}*/


//void split_by_delimiter(constconst string delim, )

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
			//TODO: Abstract out splitting string
			boost::char_separator<char> sep{","};
  			tokenizer tok{left, sep};
			vector<string> bits;
			for (const auto &b : tok)
				bits.push_back(b);
		}
		//TODO: Check bit names are valid
		// val = parse_expression(right)
	}
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
