#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include "parsing.h"
#include "Linear.h"
#include "ops.h"
#include "utils.h"

using namespace std;

map<string, Linear> varset;

struct expr {
	string l;
	string op;
	string r;
};

string clean_expr(string s) {
	boost::algorithm::trim(s);
	if (s == "" || s[0] != '(' || s.back() != ')') {
		return s;
	}
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

bool split_expr_binary(string& s, vector<string> ops, struct expr& result) {
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

Linear parse_expression(string s, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data) {
	s = clean_expr(s);
	bool split;
	if (s == "")
		throw invalid_argument("Empty expresion");
	struct expr sp;
	vector<string> delim = {"^"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l, cnts, eqs, mul_data);
		Linear right = parse_expression(sp.r, cnts, eqs, mul_data);
		Linear ret = new_xor(left, right, cnts, eqs, mul_data);
		// TODO: assert result correct
		return ret;
	}
	delim = {"+", "-"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l, cnts, eqs, mul_data);
		Linear right = parse_expression(sp.r, cnts, eqs, mul_data);
		if (sp.op == "+")
			left.add(right);
		else
			left.sub(right);
		return left;
	}
	delim = {"*", "/"};
	split = split_expr_binary(s, delim, sp);
	if (split) {
		Linear left = parse_expression(sp.l, cnts, eqs, mul_data);
		Linear right = parse_expression(sp.r, cnts, eqs, mul_data);
		if (sp.op == "*") {
			Linear ret = new_multiplication(left, right, cnts, eqs, mul_data);
			return ret;
		} else {
			Linear ret = new_division(left, right, cnts, eqs, mul_data);
			return ret;
		}
	}
	if (s.length() > 5 && s.substr(0,5) == "bool(") {
		Linear ret = parse_expression(s.substr(4), cnts, eqs, mul_data);
		Linear tmp = ret;
		Linear tmp2 = ret;
		Linear one;
		new_const(1, one);
		tmp2.sub(one);
		new_multiplication(tmp, tmp2, cnts, eqs, mul_data, false);
		cnts.bit_count += 1;
		return ret;
	}
	if (s[0] == '-') {
		Linear ret = parse_expression(s.substr(1), cnts, eqs, mul_data);
		ret.mul(mod-1);
		return ret;
	}
	if (regex_match(s, var_re)) {
		if (varset.count(s) != 0) {
			return varset[s];
		}
		else
			throw invalid_argument("Variable not defined");
	}
	if (regex_match(s, secret_re)) {
		Linear ret;
		mpz_class mpz_val = mpz_class(s.substr(1));
		new_temp(mpz_val, ret, cnts);
		return ret;
	}
	if (regex_match(s, num_re)) {
		Linear ret;
		mpz_class mpz_val = mpz_class(s);
		new_const(mpz_val, ret);
		return ret;
	}
	throw invalid_argument("Cannot parse expression");
}

vector<Linear> parse_expressions(string s, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data) {
	struct expr sp;
	vector<string> delim = {","};
	bool split = split_expr_binary(s, delim, sp);
	if (split) {
		vector<Linear> l = parse_expressions(sp.l, cnts, eqs, mul_data);
		Linear r = parse_expression(sp.r, cnts, eqs, mul_data);
		l.push_back(r);
		return l;
	}
	Linear ret_val = parse_expression(s, cnts, eqs, mul_data);
	vector<Linear> ret;
	ret.push_back(ret_val);
	return ret;
}

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
void parse_statement(vector<Linear>& eqs, string& s, struct counts& cnts, vector<struct mul>& mul_data) {
	//TODO: strip s
	if (s.length() > 6 && s.substr(0,6) == "debug "){
		Linear lin = parse_expression(s.substr(6), cnts, eqs, mul_data);
		cout << "DEBUG " << s.substr(6) << ": " << lin.real << endl;
		lin.to_str();
		cout << "const => " << lin.constant << endl;
		return;
	}
	vector<string> assn_ops = {":=", "=:", "==", "=", "=?"};
	struct expr sp;
	bool match = split_expr_binary(s, assn_ops, sp);
	if (match) {
		string left = sp.l;
		string op = sp.op;
		string right = sp.r;

		if (op == ":=") {
			vector<string> bits;
			boost::char_separator<char> sep{","};
  			tokenizer tok{left, sep};
			for (const auto &b : tok) {
				bits.push_back(b);
			}
			//TODO: Check bit names are valid
			Linear val = parse_expression(right, cnts, eqs, mul_data);
			//TODO: assert bits length <= 256
			vector<Linear> bitvars;
			if (val.is_const()) {
				for (int i = 0; i < bits.size(); i++) {
					Linear l;
					mpz_class v = (val.real >> i) & 1;
					new_const(v,l);
					bitvars.push_back(l);
				} 
			} else {
				for (int i = 0; i < bits.size(); i++) {
					Linear l;
					mpz_class v = (val.real >> i) & 1;
					new_temp(v, l, cnts);
					bitvars.push_back(l);
					Linear tmp = l;
					Linear tmp2 = l;
					Linear one;
					new_const(1, one);
					tmp2.sub(one);
					Linear eq = new_multiplication(tmp, tmp2, cnts, eqs, mul_data);
					eqs.push_back(eq);
					tmp = l;
					mpz_class shifted_i = shift_left(1, i);
					tmp.mul(shifted_i);
					val.sub(tmp);
				}
				eqs.push_back(val);
			}
			for (int i = 0; i < bits.size(); i++) {
				varset[bits[i]] = bitvars[i];
			}
		} else if (op == "=:") {
			vector<Linear> bits;
			bits = parse_expressions(right, cnts, eqs, mul_data);
			// TODO: Check varible name on left valid, len of bits <=256
			mpz_class real = 0;
			for (int i = 0; i < bits.size(); i++) {
				//TODO: Assert bit either 0 or 1
				real = real + bits[i].real * (1 << i);
				if (all_const(bits)) {
					Linear val;
					new_const(real, val);
					varset[left] = val;
				} else {
					Linear val;
					new_temp(real, val, cnts);
					varset[left] = val;
					for (int i = 0; i < bits.size(); i++) {
						Linear tmp = bits[i];
						mpz_class shifted_i = shift_left(1, i);
						tmp.mul(shifted_i);
						val.sub(tmp);
					}
					eqs.push_back(val);
				}
			}
		} else if (op == "=") {
			//TODO: check left is proper variable name
			Linear ex = parse_expression(right, cnts, eqs, mul_data);
			varset[left] = ex;
		} else if (op == "==") {
			Linear l = parse_expression(left, cnts, eqs, mul_data);
			Linear r = parse_expression(right, cnts, eqs, mul_data);
			l.sub(r);
			cout << "checking equality: " << endl;
			cout << l.real << "==" << r.real << endl;
			//TODO: check l.val == r.val
			eqs.push_back(l);
			l.to_str();
		} else if (op == "=?") {
			Linear ex = parse_expression(right, cnts, eqs, mul_data);
			int is_zero = ex.is_zero();
			Linear lin;
			lin.real = is_zero;
			lin.constant = is_zero;
			varset[left] = lin;
		}
	}
}