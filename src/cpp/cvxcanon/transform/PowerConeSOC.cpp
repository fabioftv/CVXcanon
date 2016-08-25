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
    Expression t = epi_var(expr, "power");
    Expression vec_t = reshape(t, dim(x), 1);
    Expression vec_x = reshape(x, dim(x), 1);
    constraints->push_back(
        soc(hstack({add(constant(1), neg(vec_t)), mul(constant(2), vec_x)}),
            add(constant(1), vec_t)));
    return t;
  }
  else if (p > 0 && p < 1){
    Expression t = epi_var(expr, "power");
    Expression vec_t = reshape(t, dim(x), 1);
    Expression vec_x = reshape(x, dim(x), 1);
    constraints->push_back(
        soc(hstack({add(constant(1), neg(vec_t)), mul(constant(2), vec_x)}),
            add(constant(1), vec_t)));
    return t;
  }
  else if (p == 1){
    return x;
  }
  else if (p > 1 && p != 2){
    Expression t = epi_var(expr, "power");
    Expression vec_t = reshape(t, dim(x), 1);
    Expression vec_x = reshape(x, dim(x), 1);
    constraints->push_back(
        soc(hstack({add(constant(1), neg(vec_t)), mul(constant(2), vec_x)}),
            add(constant(1), vec_t)));
    return t;
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

std::unordered_map<int, TransformFunction> kTransforms = {
  {Expression::GEO_MEAN, &transform_geo_mean},
};

Expression transform_expression(
    const Expression& expr,
    std::vector<Expression>* constraints) {
  // First transform children
  std::vector<Expression> linear_args;
  for (const Expression& arg : expr.args())
    linear_args.push_back(transform_expression(arg, constraints));

  // New expression, with linear args
  Expression output = Expression(expr.type(), linear_args, expr.attr_ptr());

  // Now transform non-linear functions, if necessary
  VLOG(2) << "transform_func: " << format_expression(expr);
  auto iter = kTransforms.find(expr.type());
  if (iter != kTransforms.end())
    output = iter->second(output, constraints);

  return output;
}

Problem LinearConeTransform::apply(const Problem& problem) {
  std::vector<Expression> constraints;
  Expression linear_objective = transform_expression(
      problem.objective, &constraints);
  for (const Expression& constr : problem.constraints) {
    constraints.push_back(
        transform_expression(constr, &constraints));
  }
  return {problem.sense, linear_objective, constraints};
}

bool LinearConeTransform::accepts(const Problem& problem) {
  for (const Expression& constr : problem.constraints)
    if (!have_transform(constr))
      return false;
  return have_transform(problem.objective);
}
