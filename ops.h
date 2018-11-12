#ifndef OPS_H
#define OPS_H

#include "Linear.h"

using namespace std;

void new_mul(mpz_class l, mpz_class r, Linear& nl, Linear& nr, Linear& no, struct counts& cnts, vector<struct mul>& mul_data);

void new_temp(mpz_class v, Linear& nt, struct counts& cnts);

void new_const(mpz_class v, Linear& nc);

Linear new_multiplication(Linear& l, Linear& r, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data, bool addeqs = true);

Linear new_division(Linear& l, Linear& r, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data);

Linear new_xor(Linear& l, Linear& r, struct counts& cnts, vector<Linear>& eqs, vector<struct mul>& mul_data);

#endif