#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include "Linear.h"


struct counts {
	unsigned int mul_count;
	unsigned int temp_count;
	unsigned int bit_count;
};

struct mul {
	mpz_class l;
	mpz_class r;
	mpz_class o;
};

const mpz_class mod = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141_mpz;
const mpz_class half_mod = (mod+1)/2;

mpz_class shift_left(mpz_class x, int i);

mpz_class modinv(mpz_class x, mpz_class mod);

unsigned long next_power_of_two(unsigned long v);

std::string var_str(int x);

int var_idx(int x);

int var_offset(char type);

char var_type(int x);

int var_idx(int x);

int random_variable(unsigned int range);

bool all_const(vector<Linear> v);

#endif