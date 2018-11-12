#ifndef PARSING_H
#define PARSING_H

#include <regex>
#include <string>
#include <vector>
#include "Linear.h"

using namespace std;

const regex var_re = regex("[A-Za-z_][0-9a-zA-Z_]*");
const regex secret_re = regex("#(-?[0-9]+)");
const regex num_re = regex("[0-9]+");

string clean_expr(string s);

bool split_expr_binary(string& s, vector<string> ops, struct expr& result);

Linear parse_expression(string s, struct counts& cnts);

vector<Linear> parse_expressions(string s, struct counts& cnts);

void parse_statement(vector<Linear>& eqs, string& s, struct counts& cnts, vector<struct mul>& mul_data);

#endif