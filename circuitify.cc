#include<stdio.h>
#include<assert.h>
#include<math.h>
#include<gmp.h>
#include<vector>
#include<tuple>
#include<map>
#include<string>
#include <stdexcept>
#include <gmpxx.h>
using namespace std;

#define MODULUS 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
#define HALFMODULUS = ((MODULUS + 1)/2)
#define COST_SCALAR_MUL 5
#define COST_SCALAR_NEG 2
#define COST_SCALAR_COPY 1

typedef tuple<char, int> var_tuple;

vector<struct mul> mul_data;
int mul_count = 0;


class Vars {
public:
	map<int, mpz_class> var_map;

/*	Vars(const Vars& obj) {
		var_map = obj.var_map;
	}*/

	void add_var(char type, int idx, char *val) {
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
		mpz_init(var_map[idx*4 + s]);
		mpz_set_str(var_map[idx*4 + s], val, 10);
	}

	void add(Vars& other) {
		//var_map.insert(other.var_map.begin(), other.var_map.end())
		/*for (auto const&x : other.var_map) {
			if (var_map.find(x.first) == var_map.end()) {
				mpz_init(var_map[x.first]);
				mpz_set(var_map[x.first], x.second);
			} else {
				mpz_add(var_map[x.first], var_map[x.first], x.second);
			}
		}*/

	}
};

class Linear {

	public:
	mpz_t val;
	bool is_const;
	Vars vars; 

	Linear(const char *v) {
		mpz_init(val);
		mpz_set_str(val, v, 10);
	}

	Linear(const mpz_t v) {
		mpz_init(val);
		mpz_set(val, v);
	}

	Linear multiply(mpz_t v) {
		mpz_t new_val;
		mpz_init(new_val);
		mpz_mul(new_val, val, v);
		Linear result(new_val);
		result.is_const = is_const;
		return result;
	}

	/*void add(Linear& other) {
		mpz_add(val, val, other.val);
		vars.add(other.vars);
	}*/
};

struct mul {
	mpz_t l;
	mpz_t r;
	mpz_t o;
};

struct multiplication {
	Linear l;
	Linear r;
	Linear o;
};

int modinv(int n) {
	assert(n != 0);
	return 5;//(int)pow(n, MODULUS-2);
}

void new_mul(mpz_t l, mpz_t r, struct linear *nl, struct linear *nr, struct linear *no)
{
  /*struct multiplication mult = {l,r,};
  mul_data.push_back(mult);
  mul_count += 1;
  struct mul m = {l.real, r.real, l.real}*/
	return;
}

Linear new_multiplication(Linear& l, Linear& r) {
	if (l.is_const)
		return r.multiply(l.val);
	if (r.is_const)
		return l.multiply(r.val);
	if (mpz_cmp(l.val, r.val) < 0) {
		Linear tmp = l;
		l = r;
		r = tmp;
	}
}


int main() {
//  printf("%d\n", modinv(5));
  
  mpz_t n;
  mpz_init(n);
  mpz_set_ui(n,0);
  printf("n = ");
  mpz_out_str(stdout, 10, n);
  printf("\n");
  return 0;
}
