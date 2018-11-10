#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <chrono>
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
#include <iostream>
#include <fstream>
#include <bitset>
#include <climits>
using namespace std;

#define COST_SCALAR_MUL 5
#define COST_SCALAR_NEG 2
#define COST_SCALAR_COPY 1

const string SECRET_FILENAME = "/Users/michaelstraka/bulletproofs_research/secp256k1-mw/src/modules/bulletproofs/bin_circuits/SHA2cpp.assn";
const string CIRCUIT_FILENAME = "/Users/michaelstraka/bulletproofs_research/secp256k1-mw/src/modules/bulletproofs/bin_circuits/SHA2cpp.circ";


typedef tuple<char, int> var_tuple;

vector<struct mul> mul_data;

unsigned int mul_count = 0;
unsigned int temp_count = 0;
unsigned int bit_count = 0;
mpz_class mod = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141_mpz;
mpz_class half_mod = (mod+1)/2;

regex var_re = regex("[A-Za-z_][0-9a-zA-Z_]*");
regex secret_re = regex("#(-?[0-9]+)");
regex num_re = regex("[0-9]+");


template <typename T>
T swap_endian(T u) {
    static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");
    union {
        T u;
        unsigned char u8[sizeof(T)];
    } source, dest;
    source.u = u;
    for (size_t k = 0; k < sizeof(T); k++)
        dest.u8[k] = source.u8[sizeof(T) - k - 1];
    return dest.u;
}

mpz_class modinv(mpz_class x) {
	mpz_t ret;
	mpz_init(ret);
	mpz_class n = mod - 2;
	mpz_powm(ret, x.get_mpz_t(), n.get_mpz_t(), mod.get_mpz_t());
	return mpz_class(ret);
}

unsigned long next_power_of_two(unsigned long v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

string var_str(int x) {
	int type = x % 4;
	int val = x / 4;
	string result;
	if (type == 0) 
		result = "L";
	else if (type == 1)
		result = "R";
	else if (type == 2)
		result = "O";
	else if (type == 3)
		result = "T";
	result.append(to_string(val));
	return result;
}

//TODO: This is bad abstraction
int random_variable() {
	int offset = rand() % 3; //'L', 'R', or 'O'
	int idx = rand() % mul_count;
	return 4*idx + offset;
}

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

char var_type(int x) {
	int offset = x % 4;
	if (offset == 0)
		return 'L';
	else if (offset == 1)
		return 'R';
	else if (offset == 2)
		return 'O';
	else if (offset == 3)
		return 'T';
	else
		throw invalid_argument("Modular arithmetic is broken");
}

int var_idx(int x) {
	return x / 4;
}

class Vars {

private:

	// deletes variables with value 0
	void reduce() {
		for(map<int, mpz_class>::iterator it = var_map.begin(); it != var_map.end();) {
			if((it->second % mod) == 0) {
				it = var_map.erase(it);
			} else {
				it++;
			}
		}
	}	

public:
	map<int, mpz_class> var_map;


	void add_var(char type, int idx, int val) {
		int s = var_offset(type);
		var_map[4*idx + s] = val % mod;
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

	void index_temp_vars(map<int, vector<int>>& index, int pos) {
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

	int num_vars() {
		return var_map.size();
	}

	void add(Vars& other) {
		for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				var_map[x.first] = x.second % mod;
			} else {
				var_map[x.first] = (var_map[x.first] + x.second) % mod;
			}
		}
		reduce();
	}

	void sub(Vars& other) {
		for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				var_map[x.first] = (mod - x.second) % mod;
			} else {
				var_map[x.first] = (var_map[x.first] - x.second + mod) % mod;
			}
		}
		reduce();
	}

	void mul(mpz_class v) {
		for (auto const&x : var_map) {
			var_map[x.first] = (var_map[x.first] * v) % mod;
		}
	}

	void div(mpz_class v) {
		for (auto const&x : var_map) {
			var_map[x.first] = (var_map[x.first] / v) % mod;
		}
	}	

	bool is_empty() {
		return var_map.empty();
	}

	bool is_zero() {
		for (auto const&x : var_map) {
			if (x.second != 0) {
				return false;
			}
		}
		return true;
	}

	/*void to_str() {
		for (std::map<int,mpz_class>::iterator it=var_map.begin(); it!=var_map.end(); ++it) {
    		print_var(it->first);
    		std::cout << " => " << it->second << '\n';
    	}
	}*/

	void str(mpz_class constant) {
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

	int cost() {
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

	map<int, mpz_class>::iterator vars_begin() {
		return var_map.begin();
	}

	map<int, mpz_class>::iterator vars_end() {
		return var_map.end();
	}
};

class Linear {

	public:
	mpz_class real;
	mpz_class constant;
	Vars vars; 

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
		real = (real + other.real) % mod;
		constant = (constant + other.constant) % mod;
		vars.add(other.vars);
	}

	void sub(Linear& other) {
		real = (real - other.real + mod) % mod;
		constant = (constant - other.constant + mod) % mod;
		vars.sub(other.vars);
	}

	void mul(mpz_class v) {
		real = (real * v) % mod;
		constant = (constant * v) % mod;
		vars.mul(v);
	}

	void div(mpz_class v) {
		v = modinv(v);
		mul(v);
	}

	void index_temp_vars(map<int, vector<int>>& index, int pos) {
		vars.index_temp_vars(index, pos);
	}

	void to_str() {
		vars.str(constant);
	}

	int equation_cost() {
		return vars.cost();
	}

	bool is_const() {
		return vars.is_empty();
	}

	bool is_zero() {
		if (constant != 0) {
			return false;
		}
		return vars.is_zero();
	}

	map<int, mpz_class>::iterator vars_begin() {
		return vars.vars_begin();
	}

	map<int, mpz_class>::iterator vars_end() {
		return vars.vars_end();
	}
};

vector<Linear> eqs;
map<string, Linear> varset;

struct mul {
	mpz_class l;
	mpz_class r;
	mpz_class o;
};

struct eq_val {
	unsigned int eq_num;
	mpz_class val;
};

struct multiplication {
	Linear l;
	Linear r;
	Linear o;
};

void new_mul(mpz_class l, mpz_class r, Linear& nl, Linear& nr, Linear& no)
{
	mpz_class o = (l*r) % mod;
	struct mul m = {l, r, o};
	mul_data.push_back(m);
	nl.real = l;
	nl.constant = 0;
	nl.add_var('L', mul_count, 1);
	nr.real = r;
	nr.constant = 0;
	nr.add_var('R', mul_count, 1);
	no.real = o;
	no.constant = 0;
	no.add_var('O', mul_count, 1);
	mul_count += 1;
}

void new_temp(mpz_class v, Linear& nt) {
	nt.real = v;
	nt.constant = 0;
	nt.add_var('T', temp_count, 1);
	temp_count += 1;
}

void new_const(mpz_class v, Linear& nc) {
	nc.real = v;
	nc.constant = v;
}

// currently mutates l and/or r
Linear new_multiplication(Linear& l, Linear& r, bool addeqs = true) {
	if (l.is_const()) {
		r.mul(l.constant);
		return r;
	}
	if (r.is_const()) {
		l.mul(r.constant);
		return l;
	}

	if (r.constant < l.constant) {
		Linear tmp = l;
		l = r;
		r = tmp;
	}
	Linear lv = Linear();
	Linear rv = Linear();
	Linear ret = Linear();
	new_mul(l.real, r.real, lv, rv, ret);
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
	if (r.is_const()) {
		l.div(r.constant);
		return l;
	}
	Linear lv = Linear();
	Linear rv = Linear();
	Linear ret = Linear();
	new_mul((l.real * modinv(r.real)) % mod, r.real, ret, rv, lv);
	l.sub(lv);
	r.sub(rv);
	eqs.push_back(l);
	eqs.push_back(r);
	return ret;
}

// mutates l and/or r
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

struct expr {
	string l;
	string op;
	string r;
};

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

Linear parse_expression(string s) {
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
		new_temp(mpz_val, ret);
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

bool all_const(vector<Linear> v) {
	for (int i = 0; i < v.size(); i++) {
		if (!v[i].is_const())
			return false;
	}
	return true;
}

mpz_class shift_left(mpz_class x, int i) {
	mpz_t result;
	mpz_init(result);
	mp_bitcnt_t bits = i;
	mpz_mul_2exp(result, x.get_mpz_t(), bits);
	return mpz_class(result);
}

typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
void parse_statement(string& s) {
	//TODO: strip s
	if (s.length() > 6 && s.substr(0,6) == "debug "){
		Linear lin = parse_expression(s.substr(6));
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
			Linear val = parse_expression(right);
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
			bits = parse_expressions(right);
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
					new_temp(real, val);
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
			Linear ex = parse_expression(right);
			varset[left] = ex;
		} else if (op == "==") {
			Linear l = parse_expression(left);
			Linear r = parse_expression(right);
			l.sub(r);
			cout << "checking equality" << endl;
			cout << l.real << " " << r.real << endl;
			//TODO: check l.val == r.val
			eqs.push_back(l);
			l.to_str();
		} else if (op == "=?") {
			Linear ex = parse_expression(right);
			int is_zero = ex.is_zero();
			Linear lin;
			lin.real = is_zero;
			lin.constant = is_zero;
			varset[left] = lin;
		}
	}
}


void pivot_variable_temp(char type, int idx, map<int, vector<int>>& index, vector<int>& to_eliminate, bool eliminate = false) {
	int c = 0;
	int cc = 0;
	int low = -1;
	Linear leq;
	vector<int> vec;
	vector<int> temp_eqs = index[idx];
//	cout << "Size: " << temp_eqs.size() << endl;
	for (int i = 0; i < temp_eqs.size(); i++) {
		/*if (find(to_eliminate.begin(), to_eliminate.end(), temp_eqs[i]) != to_eliminate.end()) {
			continue;
		}*/
		Linear &lin = eqs[temp_eqs[i]];
		vec.push_back(temp_eqs[i]);
		cc += 1;
		if (low == -1 || c > lin.num_vars()) {
			low = temp_eqs[i];
			c = lin.num_vars();
			Linear tmp = lin;
			mpz_class v = lin.get_var(type, idx);
			mpz_class inv = modinv(v);
			tmp.mul(inv);
			leq = lin;
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
	//	if (find(to_eliminate.begin(), to_eliminate.end(), low) == to_eliminate.end())
			to_eliminate.push_back(low);
	}
}

void pivot_variable(vector<Linear>& eqs_vec, int vnam, bool eliminate = false) {
	char type = var_type(vnam);
	int idx = var_idx(vnam);
	int c = 0;
	int cc = 0;
	int low = -1;
	Linear leq;
	vector<int> vec;
//	cout << "pivot_varing, type: " << type << " idx: " << idx << endl;
	for (int i = 0; i < eqs_vec.size(); i++) {
		if (eqs_vec[i].has_var(type, idx)) {
	//		cout << "found it!" << endl;
			vec.push_back(i);
			cc += 1;
			if (low == -1 || c > eqs_vec[i].num_vars()) {
				low = i;
				c = eqs_vec[i].num_vars();
				Linear tmp = eqs_vec[i];
				mpz_class v = eqs_vec[i].get_var(type, idx);
				mpz_class inv = modinv(v);
				tmp.mul(inv);
				leq = eqs_vec[i];
			}
		}
	}
//	cout << "cc: " << cc << endl;
	if (cc > 1) {
		for (auto i: vec) {
			if (i != low) {
				Linear tmp = leq;
				tmp.mul(eqs_vec[i].get_var(type, idx));
				int before = eqs_vec[i].equation_cost();
				eqs_vec[i].sub(tmp);
				int after = eqs_vec[i].equation_cost();
				int diff = after - before;
			//	cout << diff << endl;
			}
		}
	}
	if (eliminate and cc > 0) {
		eqs_vec.erase(eqs_vec.begin()+low);
	}
}

map<int, vector<int>> index_temp_vars() {
	map<int, vector<int>> temp_index;
	for (int i = 0; i < eqs.size(); i++) {
		eqs[i].index_temp_vars(temp_index, i);
	}
	return temp_index;
}

vector<Linear> eliminate_indices(vector<Linear> vec, vector<int> to_eliminate) {
	//int *a = &vec[0];
	sort(to_eliminate.begin(), to_eliminate.end());
	int j = 0;
	vector<Linear> result;
	for (int i = 0; i < vec.size(); i++) {
		if (to_eliminate[j] == i) {
			j += 1;
			continue;
		}
		result.push_back(vec[i]);
	}
	return result;
}

void eliminate_temps() {
	auto start = chrono::high_resolution_clock::now();
	map<int, vector<int>> index = index_temp_vars();
	auto finish = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Time to index temp vars: " << elapsed.count() << endl;


	
	vector<int> to_eliminate;
	for (int i = 0; i < temp_count; i++) {
		if (i % 250 == 0)
			cout << "temp_eliminated: " << i << "/" << temp_count << endl;
		pivot_variable_temp('T', i, index, to_eliminate, true);
	}
	cout << "eliminating" << endl;
	start = chrono::high_resolution_clock::now();
	sort(to_eliminate.begin(), to_eliminate.end(), greater<int>());
	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Time to sort to_eliminate: " << elapsed.count() << endl;
	start = chrono::high_resolution_clock::now();


	/*for (int x : to_eliminate) {
		eqs.erase(eqs.begin()+x);
	}*/

	eqs = eliminate_indices(eqs, to_eliminate);


	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Time to delete temp vars: " << elapsed.count() << endl;


}

void print_andytoshi_format() {
	printf("%d,0,%d,%lu; ", mul_count, bit_count, eqs.size());
	for (int i = 0; i < eqs.size(); i++) {
		int pos = 0;
		for (auto e: eqs[i].vars.var_map) {
			int key = e.first;
			mpz_class val = e.second % mod;

			bool negative = false;
			if (val*2 > mod) {
				negative = true;
				val = mod - val;
			}
			if (pos == 0 && negative)
				cout << "-";
			else if (pos > 0 && negative)
				cout << " - ";
			else if (pos > 0) 
				cout << " + ";
			if (val > 1)
				cout << val << "*";
			cout << var_str(key);
			pos += 1;
		}
		cout << " = ";
		mpz_class val = (mod - eqs[i].constant) % mod;
		if (val * 2 > mod) {
			cout << "-";
			val = mod - val;
		}
		cout << val << "; ";
	}
	cout << endl;
}

int eqs_cost(vector<Linear>& eqs_vec) {
	int cost = 0;
	for (auto& x : eqs_vec) {
		cost += x.equation_cost();
	}
	return cost;
}

void reduce_eqs(int iter) {
	auto start = chrono::high_resolution_clock::now();
	int cost = eqs_cost(eqs);
	for (int i = 0; i < iter; i++) {
		vector<Linear> neweqs = eqs;
		if (i % 10 == 0) {
			auto now = chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = now - start;
			cout << "[" << elapsed.count() << "]" << " Reduced to " << eqs_cost(eqs)
			<< " cost (step " << i << "/" << iter << ")" << endl;
		}
		for (int j = 0; j < 4; j++) {
			int vnam = random_variable();
			pivot_variable(neweqs, vnam);
			int neweqs_cost = eqs_cost(neweqs);
			if (neweqs_cost < cost) {
				eqs = neweqs;
				cost = neweqs_cost;
				break;
			}
		}
	}
	auto finish = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << "Took " << elapsed.count() << " to reduce eqs" << endl;
}

string hex_word(mpz_class val) {
	ostringstream oss;
	oss << hex << val;
	string result = oss.str();
	if (result.length() == 1) {
		result.insert(0, "0");
	}
	return result;
}

string encode_scalar_hex(mpz_class val) {
	ostringstream oss;
	if (val < 0)
		throw invalid_argument("encode_scalar_hex received negative value");

	mpz_class prefix;
	if (val > half_mod) {
		prefix = 0x80;
		val = mod - val;
	} else {
		prefix = 0x00;
	}
	int count = 0;
	while (val > 0) {
		oss << hex_word(val % 0x100);
		val = val >> 8;
		count += 1;
	}
	string result = oss.str();
	prefix = prefix | count;
	result.insert(0, hex_word(prefix));

	return result;
}

// each word is 2 bytes
string split_into_words(string hex, int n_words, bool use_padding = true) {
	string result = "";

	int padding;
	if (use_padding) {
		padding = 4 - (hex.length() % 4);
		if (padding == 4)
			padding = 0;
		padding += 4*n_words - (hex.length()+padding);
	} else {
		padding = 0;
	}
	for (int i = 0; i < padding; i++) {
		if (i != 0 && i%4 == 0)
			result += " ";
		result += "0";
	}
	for (int i = 0; i < hex.length(); i++) {
		if (i != 0 && (i+padding) % 4 == 0)
			result += " ";
		result += hex[i];
	}
	return result;
}

template <typename T>
string hex_string(T val) {
	ostringstream oss;
	oss << hex << val;
	return oss.str();
}

template <typename T>
void write_enc(T val, ofstream& f) {
	T i = swap_endian<T>(val);
	string h = hex_string(i);
	f << split_into_words(h, sizeof(T)/2) << " ";
}

// does not pad
void write_mpz_enc(mpz_class val, ofstream& f) {
	f << split_into_words(encode_scalar_hex(val), 0, false) << " ";
}

void write_secret_data(string filename) {
	ofstream f;
	f.open(filename);
	// 2 bytes version (1), 2 bytes flags (0), 4 bytes n_commits (0), 8 bytes n_gates
	write_enc<uint32_t>(1, f);
	write_enc<uint32_t>(0, f);
	write_enc<uint64_t>(mul_count, f);

	for (int i = 0; i < mul_data.size(); i++) {
		write_mpz_enc(mul_data[i].l, f);
	}
	for (int i = 0; i < mul_data.size(); i++) {
		write_mpz_enc(mul_data[i].r, f);
	}
	for (int i = 0; i < mul_data.size(); i++) {
		write_mpz_enc(mul_data[i].o, f);
	}

	f.close();
}

void write_matrix(map<int, vector<struct eq_val>>& mat, ofstream& f) {
	for (auto& e : mat) {
		vector<struct eq_val> col = e.second;
		write_enc<uint32_t>(col.size(), f);
		for (int i = 0; i < col.size(); i++) {
			write_enc<uint32_t>(col[i].eq_num, f);
			write_mpz_enc(col[i].val, f);
		}
	}
}

void write_circuit_data(string filename) {
	ofstream f;
	f.open(filename);
	unsigned int next_mul_count = next_power_of_two(mul_count);
	// 2 bytes version (1), 2 bytes flags (0), 4 bytes n_commits (0), 8 bytes n_gates, 8 bytes n_bits, 8 bytes n_constraints
	write_enc<uint32_t>(1, f);
	write_enc<uint32_t>(0, f);
	write_enc<uint64_t>(next_mul_count, f);
	write_enc<uint64_t>(bit_count, f);
	write_enc<uint64_t>(eqs.size(), f);

	map<int, vector<struct eq_val>> WL;
	map<int, vector<struct eq_val>> WR;
	map<int, vector<struct eq_val>> WO;
	vector<mpz_class> C;
	for (unsigned int i = 0; i < eqs.size(); i++) {
		for (auto it = eqs[i].vars_begin(); it != eqs[i].vars_end(); ++it) {
			char type = var_type(it->first);
			int idx = var_idx(it->first);

			struct eq_val ev = {i, it->second};
			if (type == 'L') {
				WL[idx].push_back(ev);
			} else if (type == 'R') {
				WR[idx].push_back(ev);
			} else if (type == 'O') {
				WO[idx].push_back(ev);
			} else {
				throw invalid_argument("unknown type encountered in writing circuit data");
			}
		}
	/*	for (int i = 0; i < (next_mul_count-mul_count); i++) {
			write_enc<uint32_t>(0, f);
		}*/
		C.push_back(eqs[i].constant);
	}

	write_matrix(WL, f);
	write_matrix(WR, f);
	write_matrix(WO, f);
	for (int i = 0; i < C.size(); i++) {
		write_mpz_enc((mod-C[i]) % mod, f);
	}

	f.close();
}

int main() {

	string input_line;
	auto start = chrono::high_resolution_clock::now();
	while (getline(cin, input_line)) {
		cout << "input line: " << input_line << endl;
		parse_statement(input_line);
	}
	auto finish = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << endl << "Time to parse input: " << elapsed.count() << endl;


	//print_andytoshi_format();
	
	printf("%d multiplications, %d temporaries, %lu constraints, %d cost\n", mul_count, temp_count, eqs.size(), eqs_cost(eqs));

	start = chrono::high_resolution_clock::now();
	eliminate_temps();
	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Time to eliminate vars: " << elapsed.count() << endl;


	//print_andytoshi_format();

/*	start = chrono::high_resolution_clock::now();
	reduce_eqs(mul_count);
	finish = chrono::high_resolution_clock::now();
	printf("%d multiplications, %d temporaries, %lu constraints, %d cost\n", mul_count, temp_count, eqs.size(), eqs_cost(eqs));
*/
//	print_andytoshi_format();



	/*for (int i = 0; i < eqs.size(); i++) {
		for (auto it = eps[i].vars_begin(); it != eqs[i].vars_end(); ++it) {
			char type = var_type(it.first);
			int idx = var_idx(it.first);

		}
		// for index, value in eqs[i].vars:
		// W[index].push_back(i, value)
		//C.push_back(eqs[i].constant)
	}*/

	write_secret_data(SECRET_FILENAME);
	write_circuit_data(CIRCUIT_FILENAME);

	/*ofstream myfile;
	myfile.open("./testfile");
	uint32_t i = swap_endian<uint32_t>(52);
	cout << i << endl;
	//string bin = split_into_bytes(bitset<8>(i).to_string());
	cout << hex << i << endl;
	//myfile.write(bin);
	//myfile << bin;
	myfile << hex << i << endl;*/
	//cout << endl << "I wrote to a file!" << endl;
	//myfile.close();
 	return 0;
}
