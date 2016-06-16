// Utility functions for constructing and manipulating expression trees.

#ifndef CVXCANON_EXPRESSION_EXPRESSION_UTIL_H
#define CVXCANON_EXPRESSION_EXPRESSION_UTIL_H

#include <string>
#include <vector>

#include "cvxcanon/expression/Expression.hpp"

// Syntatic sugar that allows for writing expresison trees in natural functional
// form, for example:
// Expression expr = add(mul(const(A), var(10, 1)), neg(const(b)));
Expression abs(Expression x);
Expression add(Expression x, Expression y);
Expression constant(DenseMatrix value);
Expression constant(double value);
Expression hstack(std::vector<Expression> args);
Expression mul(Expression x, Expression y);
Expression neg(Expression x);
Expression p_norm(Expression x, double p);
Expression power(Expression x, double p);
Expression quad_over_lin(Expression x, Expression y);
Expression reshape(Expression x, int m, int n);
Expression sum_entries(Expression x);
Expression var(int m, int n);
Expression vstack(std::vector<Expression> args);

// Syntatic sugar for constraints
Expression eq(Expression x, Expression y);  // x == y
Expression leq(Expression x, Expression y);  // x <= y
Expression soc(Expression x, Expression y);  // ||x||_2 <= y
// K = {(x,y,z) | y > 0, ye^(x/y) <= z} U {(x,y,z) | x <= 0, y = 0, z >= 0}
Expression exp_cone(Expression x, Expression y, Expression z);

// Construct a new VAR expression of the same size as f, intended to be used in
// canonicalization when adding epigraph constraints.
Expression epi_var(const Expression& f, const std::string& name);

// Returns true if size represents a scalar, i.e. all dimensinos are 1
bool is_scalar(const Size& size);

#endif  // CVXCANON_EXPRESSION_EXPRESSION_UTIL_H
