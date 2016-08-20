#include "PowerConeTransform.hpp"

#include <unordered_map>
#include <vector>

#include "cvxcanon/expression/Expression.hpp"
#include "cvxcanon/expression/ExpressionShape.hpp"
#include "cvxcanon/expression/ExpressionUtil.hpp"
#include "cvxcanon/expression/TextFormat.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include "glog/logging.h"

typedef Expression(*TransformFunction)(
    const Expression& expr,
    std::vector<Expression>* constraints);

Expression transform_power(
    const Expression& expr,
    std::vector<Expression>* constraints) {
  const double p = expr.attr<PowerAttributes>().p;
  const Expression& x = expr.arg(0);
  if (p < 0){

  }
  else if (p > 0 && p < 1){

  }
  else if (p == 1){
    return x;
  }
  else if (p > 1 && p != 2){

  }
  else if (p == 2) {
    Expression t = epi_var(expr, "power_2");
    Expression vec_t = reshape(t, dim(x), 1);
    Expression vec_x = reshape(x, dim(x), 1);
    constraints->push_back(
        soc(hstack({add(constant(1), neg(vec_t)), mul(constant(2), vec_x)}),
            add(constant(1), vec_t)));
    return t;
  }
}
