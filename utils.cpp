#include <string>
#include <stdexcept> 
#include <stdlib.h>
#include <gmpxx.h>
using namespace std;

mpz_class modinv(mpz_class x, mpz_class mod) {
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

std::string var_str(int x) {
	int type = x % 4;
	int val = x / 4;
	std::string result;
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

int var_idx(int x) {
	return x / 4;
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

//TODO: This is bad abstraction
int random_variable(unsigned int range) {
	int offset = rand() % 3; //'L', 'R', or 'O'
	int idx = rand() % range;
	return 4*idx + offset;
}