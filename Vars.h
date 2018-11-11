#ifndef VARS_H
#define VARS_H

#include <map>
#include <gmpxx.h>

class Vars {

private:

	// deletes variables with value 0
	void reduce();	

public:
	std::map<int, mpz_class> var_map;
	mpz_class mod = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141_mpz;

	void add_var(char type, int idx, int val);

	bool has_var(char type, int idx);

	mpz_class get_var(char type, int idx);

	void index_temp_vars(std::map<int, std::vector<int>>& index, int pos);

	int num_vars();

	void add(Vars& other);

	void sub(Vars& other);

	void mul(mpz_class v);

	void div(mpz_class v);	

	bool is_empty();

	bool is_zero();

	void str(mpz_class constant);

	int cost();

	std::map<int, mpz_class>::iterator vars_begin();

	std::map<int, mpz_class>::iterator vars_end();
};

#endif