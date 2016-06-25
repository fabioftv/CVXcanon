
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



//  EXPECT_EQ("var", format_expression(x));
}

TEST(TextFormatTest2, Second) {
  Expression x = constant(10.0);
  EXPECT_EQ("const", format_expression(x));
}
