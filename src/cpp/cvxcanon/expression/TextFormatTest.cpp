
#include <string>

#include "gtest/gtest.h"

#include "cvxcanon/expression/ExpressionUtil.hpp"
#include "cvxcanon/expression/TextFormat.hpp"

extern std::unordered_map<int, std::string> kExpressionNames;
TEST(TextFormatTest, Names) {
  for (int i = 0; i < Expression::NUM_TYPES; i++) {
    EXPECT_TRUE(kExpressionNames.find(i) != kExpressionNames.end());
  }
}

TEST(TextFormatTest, FormatExpression) {
  Expression x = var(10, 5, 0);
  Expression y = var(20, 10, 0);
  Expression z = var(30, 15, 0);
  double d = 10.0;

  Expression add_vars = add(x, y);
  EXPECT_EQ("add(var, var)", format_expression(add_vars));
 
  Expression mult_vars = mul(x, y);
  EXPECT_EQ("mul(var, var)", format_expression(mult_vars));

  Expression neg_var = neg(x);
  EXPECT_EQ("neg(var)", format_expression(neg_var));

  Expression abs_var = abs(x);
  EXPECT_EQ("abs(var)", format_expression(abs_var));

  Expression trans_var = transpose(x);
  EXPECT_EQ("transpose(var)", format_expression(trans_var));

  Expression log_var = log(x);
  EXPECT_EQ("log(var)", format_expression(log_var));

  Expression exp_var = exp(x);
  EXPECT_EQ("exp(var)", format_expression(exp_var));

  Expression upper_var = upper_tri(x);
  EXPECT_EQ("upper_tri(var)", format_expression(upper_var));

  Expression diag_var = diag_vec(x);
  EXPECT_EQ("diag_vec(var)", format_expression(diag_var));

  Expression diagm_var = diag_mat(x);
  EXPECT_EQ("diag_mat(var)", format_expression(diagm_var));

  Expression trace_var = trace(x);
  EXPECT_EQ("trace(var)", format_expression(trace_var));

  Expression quad_vars = quad_over_lin(x, y);
  EXPECT_EQ("quad_over_lin(var, var)", format_expression(quad_vars));

  Expression pnorm_vars = p_norm(x, d);
  EXPECT_EQ("p_norm(var)", format_expression(pnorm_vars));

  Expression eq_vars = eq(x, y);
  EXPECT_EQ("eq(var, var)", format_expression(eq_vars));

  Expression leq_vars = leq(x, y);
  EXPECT_EQ("leq(var, var)", format_expression(leq_vars));

  Expression soc_vars = soc(x, y);
  EXPECT_EQ("soc(var, var)", format_expression(soc_vars));

  Expression expc_vars = exp_cone(x, y, z);
  EXPECT_EQ("exp_cone(var, var, var)", format_expression(expc_vars));

  Expression sdp_var = sdp(x);
  EXPECT_EQ("sdp(var)", format_expression(sdp_var));






/*

Expression power(Expression x, double p) {
  auto attr = std::make_shared<PowerAttributes>();
  attr->p = p;
  return {Expression::POWER, {x}, attr};
}

Expression reshape(Expression x, int m, int n) {
  auto attr = std::make_shared<ReshapeAttributes>();
  attr->size = {{m, n}};
  return {Expression::RESHAPE, {x}, attr};
}

Expression sum_entries(Expression x, int axis) {
  auto attr = std::make_shared<SumEntriesAttributes>();
  attr->axis = axis;
  return {Expression::SUM_ENTRIES, {x}, attr};
}

Expression index(
    Expression x, int start_i, int stop_i, int start_j, int stop_j) {
  auto attr = std::make_shared<IndexAttributes>();
  attr->keys.push_back({start_i, stop_i, 1});
  attr->keys.push_back({start_j, stop_j, 1});
  return {Expression::INDEX, {x}, attr};
}

Expression var(int m, int n, int var_id) {
  auto attr = std::make_shared<VarAttributes>();
  attr->id = var_id;
  attr->size = {{m, n}};
  return {Expression::VAR, {}, attr};
}

Expression constant(double value) {
  return constant(DenseMatrix::Constant(1, 1, value));
}

Expression constant(DenseMatrix value) {
  auto attr = std::make_shared<ConstAttributes>();
  attr->constant.dense_data = value;
  attr->constant.sparse = false;
  return {Expression::CONST, {}, attr};
}



Expression epi_var(const Expression& x, const std::string& name) {
  return epi_var_size(x, name, size(x));
}

Expression epi_var_size(
    const Expression& x, const std::string& name, Size size) {
  int var_id = rand();  // NOLINT(runtime/threadsafe_fn)
  VLOG(2) << "epi_var " << var_id << ", "
          << size.dims[0] << " x " << size.dims[1];
  return var(size.dims[0], size.dims[1], var_id);
}

bool is_scalar(const Size& size) {
  return size.dims[0] == 1 && size.dims[1] == 1;
}

*/

/*


Expression hstack(std::vector<Expression> args) {
  return {Expression::HSTACK, args};
}

Expression vstack(std::vector<Expression> args) {
  return {Expression::VSTACK, args};
}
*/
/*
  // Linear functions
  {Expression::HSTACK, "hstack"},
  {Expression::INDEX, "index"},
  {Expression::KRON, "kron"},
  {Expression::RESHAPE, "reshape"},
  {Expression::SUM_ENTRIES, "sum_entries"},
  {Expression::VSTACK, "vstack"},

  // Elementwise functions
  {Expression::ENTR, "entr"},
  {Expression::HUBER, "huber"},
  {Expression::KL_DIV, "kl_div"},
  {Expression::LOG1P, "log1p"},
  {Expression::LOGISTIC, "logistic"},
  {Expression::MAX_ELEMWISE, "max_elemwise"},
  {Expression::POWER, "power"},

  // General nonlinear functions
  {Expression::GEO_MEAN, "geo_mean"},
  {Expression::LAMBDA_MAX, "lambda_max"},
  {Expression::LOG_DET, "log_det"},
  {Expression::LOG_SUM_EXP, "log_sum_exp"},
  {Expression::MATRIX_FRAC, "matrix_frac"},
  {Expression::MAX_ENTRIES, "max_entries"},
  {Expression::NORM_NUC, "norm_nuc"},
  {Expression::SIGMA_MAX, "sigma_max"},
  {Expression::SUM_LARGEST, "sum_largest"},

  // Constraints

  {Expression::SDP_VEC, "sdp_vec"},

  // Leaf nodes
  {Expression::CONST, "const"},
  {Expression::PARAM, "param"},
  {Expression::VAR, "var"},

*/


//  EXPECT_EQ("var", format_expression(x));
}

TEST(TextFormatTest2, Second) {
  Expression x = constant(10.0);
  EXPECT_EQ("const", format_expression(x));
}
