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
Expression diag_vec(Expression x);
Expression diag_mat(Expression A);
Expression exp(Expression x);
Expression hstack(std::vector<Expression> args);
Expression index(
    Expression x, int start_i, int stop_i, int start_j, int stop_j);
Expression log(Expression x);
Expression mul(Expression x, Expression y);
Expression neg(Expression x);
Expression p_norm(Expression x, double p);
Expression power(Expression x, double p);
Expression quad_over_lin(Expression x, Expression y);
Expression reshape(Expression x, int m, int n);
Expression sum_entries(Expression x, int axis = kNoAxis);
Expression trace(Expression x);
Expression transpose(Expression x);
Expression upper_tri(Expression x);
Expression var(int m, int n);
Expression vstack(std::vector<Expression> args);

// Return expressions that promote a scalar for use in either matrix addition or
// multiplication. If no promotion is necessary, than these return the input.
Expression promote_ones(Expression x, const Size& size);  // alpha*11^T
Expression promote_identity(Expression x, int n);         // alpha*I

// Syntatic sugar for constraints
Expression eq(Expression x, Expression y);  // x == y
Expression leq(Expression x, Expression y);  // x <= y

// ||x||_2 <= y
//
// If y represents a vector, than x should be a matrix so that the constraint
// can be applied to each (x_i, y_i) where x_i is a row of x.
Expression soc(Expression x, Expression y);

// K = {(x,y,z) | y > 0, ye^(x/y) <= z} U {(x,y,z) | x <= 0, y = 0, z >= 0}
Expression exp_cone(Expression x, Expression y, Expression z);

// X >> 0
Expression sdp(Expression X);

// Construct a new VAR expression of the same size as f or specified size; used
// in canonicalization when adding epigraph constraints.
Expression epi_var(const Expression& f, const std::string& name);
Expression epi_var_size(
    const Expression& f, const std::string& name, Size size);

// Returns true if size represents a scalar, i.e. all dimensinos are 1
bool is_scalar(const Size& size);

// Basic classes of expressions, note that these are functions only of the expr
// itself, not its children.
bool is_constraint(const Expression& expr);
bool is_linear(const Expression& expr);
bool is_leaf(const Expression& expr);


#endif  // CVXCANON_EXPRESSION_EXPRESSION_UTIL_H
