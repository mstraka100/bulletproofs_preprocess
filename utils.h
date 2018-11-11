#ifndef UTILS_H
#define UTILS_H

mpz_class modinv(mpz_class x, mpz_class mod);

unsigned long next_power_of_two(unsigned long v);

std::string var_str(int x);

int var_idx(int x);

int var_offset(char type);

char var_type(int x);

int var_idx(int x);

int random_variable(unsigned int range);

#endif