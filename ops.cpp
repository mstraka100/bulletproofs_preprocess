#include "ops.h"
#include "utils.h"

void new_mul(mpz_class l, mpz_class r, Linear& nl, Linear& nr, Linear& no, struct counts& cnts, vector<struct mul>& mul_data) {
	mpz_class o = (l*r) % mod;
	struct mul m = {l, r, o};
	mul_data.push_back(m);
	nl.real = l;
	nl.constant = 0;
	nl.add_var('L', cnts.mul_count, 1);
	nr.real = r;
	nr.constant = 0;
	nr.add_var('R', cnts.mul_count, 1);
	no.real = o;
	no.constant = 0;
	no.add_var('O', cnts.mul_count, 1);
	cnts.mul_count += 1;
}

void new_temp(mpz_class v, Linear& nt, struct counts& cnts) {
	nt.real = v;
	nt.constant = 0;
	nt.add_var('T', cnts.temp_count, 1);
	cnts.temp_count += 1;
}

void new_const(mpz_class v, Linear& nc) {
	nc.real = v;
	nc.constant = v;
}

// mutates l and/or r
Linear new_multiplication(Linear& l, Linear& r, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data, bool addeqs) {
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
	new_mul(l.real, r.real, lv, rv, ret, cnts, mul_data);
	l.sub(lv);
	eqs.push_back(l);
	if (addeqs){
		r.sub(rv);
		eqs.push_back(r);
	}
	return ret;
}

// mutates l and/or r
Linear new_division(Linear& l, Linear& r, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data) {
	if (r.is_const()) {
		l.div(r.constant);
		return l;
	}
	Linear lv = Linear();
	Linear rv = Linear();
	Linear ret = Linear();
	new_mul((l.real * modinv(r.real, mod)) % mod, r.real, ret, rv, lv, cnts, mul_data);
	l.sub(lv);
	r.sub(rv);
	eqs.push_back(l);
	eqs.push_back(r);
	return ret;
}

// mutates l and/or r
Linear new_xor(Linear& l, Linear& r, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data) {
	Linear lv = Linear();
	Linear rv = Linear();
	Linear mul = Linear();
	new_multiplication(l, r, cnts, eqs, mul_data);
	l.add(r);
	mul.mul(2);
	l.sub(mul);
	return l;
}