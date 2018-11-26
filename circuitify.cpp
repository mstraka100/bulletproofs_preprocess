#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <gmpxx.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <stdio.h>
#include <assert.h>
#include "utils.h"
#include "Vars.h"
#include "Linear.h"
#include "parsing.h"
using namespace std;

const string SECRET_FILENAME = "/CircuitOutput/filename.assn";
const string CIRCUIT_FILENAME = "/CircuitOutput/filename.circ";

/* Coefficient of some variable with index of the equation it appears in */
struct eq_val {
	unsigned int eq_num;
	mpz_class val;
};

/* Given int idx, eliminates the variable of the form T<idx> from eqs. Marks for future removal by adding to to_eliminate. */
void pivot_variable_temp(vector<Linear>& eqs, int idx, map<int, vector<int>>& index, bool eliminate = false) {
	if (index.count(idx) == 0)
		return;

	int c = 0;
	int cc = 0;
	int low = -1;
	int low_idx;
	int leq;
	Linear leq_lin;
	vector<int> vec;

	cout << "pivoting idx: " << idx << endl;

	vector<int>& temp_eqs = index[idx];
	for (int i = 0; i < temp_eqs.size(); i++) {
		int lin = temp_eqs[i];
		cc += 1;
		if (eqs[lin].eliminated)
			continue;
		if (!eqs[lin].has_var('T', idx))
			continue;
		if (low == -1 || c > eqs[lin].num_vars()) {
			low_idx = i;
			c = eqs[lin].num_vars();
			leq = lin;
		}
	}

	if (cc > 1) {
		leq_lin = eqs[leq];
		mpz_class v = eqs[leq].get_var('T', idx);
		mpz_class inv = modinv(v, mod);
		leq_lin.mul(inv);
		// TODO: Parallelize this loop
		for (auto i: temp_eqs) {
			if (i != leq) {
				if (!eqs[i].has_var('T', idx)) {
					continue;
				}
				eqs[i].assign_temp_vars(leq_lin, index, i);
				Linear tmp = leq_lin;
				tmp.mul(eqs[i].get_var('T', idx));
				eqs[i].sub(tmp);
			}
		}
	}
	if (eliminate and cc > 0) {
		eqs[leq].eliminated = true;
	}
	index.erase(idx);
}

/*void pivot_variable(vector<Linear>& eqs_vec, int vnam, bool eliminate = false) {
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
				mpz_class inv = modinv(v, mod);
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
}*/

map<int, vector<int>> index_temp_vars(vector<Linear>& eqs) {
	map<int, vector<int>> temp_index;
	for (int i = 0; i < eqs.size(); i++) {
		eqs[i].index_temp_vars(temp_index, i);
	}
	return temp_index;
}

void delete_eliminated_eqs(vector<Linear>& eqs) {
	for (int i = eqs.size()-1; i >= 0; i--) {
		if (eqs[i].eliminated) {
			eqs.erase(eqs.begin() + i);
		}
	}
}

void print_andytoshi_format(vector<Linear>& eqs, struct counts& cnts) {
	printf("%d,0,%d,%lu; ", cnts.mul_count, cnts.bit_count, eqs.size());
	for (int i = 0; i < eqs.size(); i++) {
		if (eqs[i].eliminated)
			continue;
	//	cout << i << ": ";
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

void eliminate_temps(vector<Linear>& eqs, struct counts& cnts) {
	map<int, vector<int>> index = index_temp_vars(eqs);
	for (int i = 0; i < cnts.temp_count; i++) {
		if (i % 250 == 0)
			cout << "temp_eliminated: " << i << "/" << cnts.temp_count << endl;
		pivot_variable_temp(eqs, i, index, true);
	}
	delete_eliminated_eqs(eqs);
}



int eqs_cost(vector<Linear>& eqs_vec) {
	int cost = 0;
	for (auto& x : eqs_vec) {
		cost += x.equation_cost();
	}
	return cost;
}

/*void reduce_eqs(int iter, struct counts& cnts) {
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
			int vnam = random_variable(cnts.mul_count);
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
}*/

string hex_word(mpz_class val) {
	ostringstream oss;
	oss << hex << val;
	string result = oss.str();
	if (result.length() == 1) {
		result.insert(0, "0");
	}
	return result;
}

vector<char> encode_scalar_to_hex(mpz_class val) {
	vector<char> result;
	int neg_count = 0;
	while (val < 0) {
		cout << "encountered negative value when writing data: " << neg_count << endl;
		val = mod + val;
		neg_count += 1;
	}

	size_t prefix;
	if (val > half_mod) {
		prefix = 0x80;
		val = mod - val;
	} else {
		prefix = 0x00;
	}
	int count = 0;
	while (val > 0) {
		mpz_class x = val % 0x100;
		result.push_back(x.get_ui());
		val = val >> 8;
		count += 1;
	}
	prefix = prefix | count;
	result.insert(result.begin(), prefix);

	return result;
}

void write_mpz_enc(mpz_class val, ofstream& f) {
	vector<char> vec = encode_scalar_to_hex(val);
	f.write(&vec[0], vec.size());
}

template <typename T>
void write_enc_to(T val, size_t width, ofstream& f) {
	char buf[width];
	size_t i;
    for (i = 0; i < width; i++) {
         buf[i] = val;
         val >>= 8;
    }
	f.write(buf, width);
}


void write_secret_data(string filename, vector<struct mul>& mul_data, struct counts& cnts) {
	ofstream f;
	f.open(filename);

	// 2 bytes version (1), 2 bytes flags (0), 4 bytes n_commits (0), 8 bytes n_gates
	write_enc_to<uint16_t>(1, 4, f);
	write_enc_to<uint16_t>(0, 4, f);
	write_enc_to<uint64_t>(cnts.mul_count, 8, f);

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

void write_matrix(map<int, vector<struct eq_val>>& mat, size_t metadata_len, ofstream& f) {
	for (auto& e : mat) {
		vector<struct eq_val> col = e.second;
		write_enc_to<uint64_t>(col.size(), metadata_len, f);
		for (int i = 0; i < col.size(); i++) {
			write_enc_to<uint64_t>(col[i].eq_num, metadata_len, f);
			write_mpz_enc(col[i].val, f);
		}
	}
}

// returns number of bytes matrix value metadata should be encoded as
size_t encoding_length(size_t n_elems) {
	if (n_elems <= 255)
		return 1;
	else if (n_elems <= (32767 * 2 + 1))
		return 2;
	else
		return 4;
}

void write_circuit_data(string filename, vector<Linear>& eqs, struct counts& cnts) {
	ofstream f;
	f.open(filename);
	size_t metadata_len = encoding_length(eqs.size());
	unsigned int next_mul_count = next_power_of_two(cnts.mul_count);

	// 2 bytes version (1), 2 bytes flags (0), 4 bytes n_commits (0), 8 bytes n_gates, 8 bytes n_bits, 8 bytes n_constraints
	write_enc_to<uint16_t>(1, 4, f);
	write_enc_to<uint16_t>(0, 4, f);

	write_enc_to<uint64_t>(next_mul_count, 8, f);
	write_enc_to<uint64_t>(cnts.bit_count, 8, f);
	write_enc_to<uint64_t>(eqs.size(), 8, f);

	map<int, vector<struct eq_val>> WL;
	map<int, vector<struct eq_val>> WR;
	map<int, vector<struct eq_val>> WO;

	for (size_t i = 0; i < cnts.mul_count; i++) {
		WL[i] = vector<struct eq_val>();
		WR[i] = vector<struct eq_val>();
		WO[i] = vector<struct eq_val>();
	}
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
		for (int i = 0; i < (next_mul_count-cnts.mul_count); i++) {
			write_enc_to<uint32_t>(0, metadata_len, f);
		}
		C.push_back(eqs[i].constant);
	}

	write_matrix(WL, metadata_len, f);
	write_matrix(WR, metadata_len, f);
	write_matrix(WO, metadata_len, f);
	for (int i = 0; i < C.size(); i++) {
		write_mpz_enc((mod-C[i]) % mod, f);
	}

	f.close();
}

int main() {

	struct counts cnts = {0, 0, 0};

	vector<Linear> eqs;
	vector<struct mul> mul_data;
	string input_line;
	auto start = chrono::high_resolution_clock::now();
	while (getline(cin, input_line)) {
		cout << "input line: " << input_line << endl;
		parse_statement(eqs, input_line, cnts, mul_data);
	}
	auto finish = chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << endl << "Time to parse input: " << elapsed.count() << endl;

	print_andytoshi_format(eqs, cnts);
	
	printf("%d multiplications, %d temporaries, %lu constraints, %d cost\n", cnts.mul_count, cnts.temp_count, eqs.size(), eqs_cost(eqs));

	start = chrono::high_resolution_clock::now();
	eliminate_temps(eqs, cnts);
	finish = chrono::high_resolution_clock::now();
	elapsed = finish - start;
	cout << "Time to eliminate vars: " << elapsed.count() << endl;

	print_andytoshi_format(eqs, cnts);

/*	start = chrono::high_resolution_clock::now();
	reduce_eqs(mul_count);
	finish = chrono::high_resolution_clock::now();
	printf("%d multiplications, %d temporaries, %lu constraints, %d cost\n", mul_count, temp_count, eqs.size(), eqs_cost(eqs));
	print_andytoshi_format(eqs, cnts); */

	write_secret_data(SECRET_FILENAME, mul_data, cnts);
	write_circuit_data(CIRCUIT_FILENAME, eqs, cnts);

 	return 0;
}
