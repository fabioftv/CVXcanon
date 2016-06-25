
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "cvxcanon/expression/ExpressionUtil.hpp"
#include "cvxcanon/expression/TextFormat.hpp"

extern std::unordered_map<int, std::string> kExpressionNames;

TEST(TextFormatTest, Sense) {
  for (int i = 0; i <= Problem::MINIMIZE; i++) {
     EXPECT_TRUE(kSenseNames.find(i) != kSenseNames.end());
  }
}

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
  DenseMatrix A;
  std::vector<Expression> v;

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

  Expression const_var = constant(A);
  EXPECT_EQ("const", format_expression(const_var));

  Expression hstack_var = hstack(v);
  EXPECT_EQ("hstack", format_expression(hstack_var));

  Expression vstack_var = vstack(v);
  EXPECT_EQ("vstack", format_expression(vstack_var));

//TODO (fabioftv): Missing "kron", "entr", "huber", "kl_div", "log1p", "logistic", "max_elemwise", "geo_mean", "lambda_max", "log_det", "log_sum_exp", "matrix_frac", "max_entries", "norm_nuc", "sigma_max", "sum_largest", "sdp_vec", "param"
}
