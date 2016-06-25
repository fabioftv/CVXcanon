
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

  Expression ad = add(x, y);
  Expression mu = mul(x, y);
  Expression ne = neg(x);

  EXPECT_EQ("add(var, var)", format_expression(ad));
  EXPECT_EQ("mul(var, var)", format_expression(mu));
  EXPECT_EQ("neg", format_expression(ne));


//  EXPECT_EQ("var", format_expression(x));
}

TEST(TextFormatTest2, Second) {
  Expression x = constant(10.0);
  EXPECT_EQ("const", format_expression(x));
}
