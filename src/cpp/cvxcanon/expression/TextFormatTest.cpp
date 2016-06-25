
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
  int m = 10;
  int n = 5;
  int i = 15;
  int j = 20;
  std::string t = "test";

  EXPECT_EQ("var", format_expression(x));

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

  Expression reshape_vars = reshape(x, m, n);
  EXPECT_EQ("reshape(var)", format_expression(reshape_vars));

  Expression power_vars = power(x, d);
  EXPECT_EQ("power(var)", format_expression(power_vars));

  Expression sum_vars = sum_entries(x, m);
  EXPECT_EQ("sum_entries(var)", format_expression(sum_vars));

  Expression index_vars = index(x, m, n, i, j);
  EXPECT_EQ("index(var)", format_expression(index_vars));

  Expression epi = epi_var(x, t);
  EXPECT_EQ("epi_var_size(var, test, 10)", format_expression(epi));

/*

Expression epi_var(const Expression& x, const std::string& name) {
  return epi_var_size(x, name, size(x));
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
  {Expression::KRON, "kron"},
  {Expression::VSTACK, "vstack"},

  // Elementwise functions
  {Expression::ENTR, "entr"},
  {Expression::HUBER, "huber"},
  {Expression::KL_DIV, "kl_div"},
  {Expression::LOG1P, "log1p"},
  {Expression::LOGISTIC, "logistic"},
  {Expression::MAX_ELEMWISE, "max_elemwise"},


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


*/



}

TEST(TextFormatTest2, Second) {
  Expression x = constant(10.0);
  EXPECT_EQ("const", format_expression(x));
}
