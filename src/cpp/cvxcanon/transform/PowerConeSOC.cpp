#include "PowerConeSOC.hpp"

#include <unordered_map>
#include <vector>

#include "cvxcanon/transform/PowerGeoMean.hpp"
#include "cvxcanon/expression/Expression.hpp"
#include "cvxcanon/expression/ExpressionShape.hpp"
#include "cvxcanon/expression/ExpressionUtil.hpp"
#include "cvxcanon/expression/TextFormat.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include "glog/logging.h"

PowerConeSOCTransform::PowerConeSOCTransform() {}
PowerConeSOCTransform::~PowerConeSOCTransform() {}

std::vector<Expression> PowerConeSOCTransform::gm_constrs(const Expression& t, 
   const Expression& expr, std::vector<std::pair<double, double>> p) {

   GeoMeanIneq geo_mean;

   assert (geo_mean.weight_vector_test(p) == true);

   std::vector<std::pair<double, double>> w(p.size());
   w = geo_mean.dyad_completion(p);

   std::vector<std::pair<std::vector<std::pair<double, double>>, 
      std::pair<std::vector<std::pair<double, double>>, 
         std::vector<std::pair<double, double>>>>> tree;

   tree = geo_mean.decompose(w);

//   t.attr<VarAttributes>().size;


/*
   std::vector<double> d;
   d(w) = t;

   if (variables.size() < w.size()){
      variables += t;
   }
   variables.size = w.size;
   
   std::vector<int> tmp;

   for (int i = 0; i < p.size; i++){
      if (p[i] > 0){
         if (i = 0){
            tmp[i] = size.w;
         }
         else{
            tmp[i] = 1;
         }
      }
      d[tpm[i]] = v;
   }
*/

   std::vector<Expression> constraints;
   return constraints;
}


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
  {Expression::GEO_MEAN_INEQ, &transform_geo_mean_ineq},
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

Problem PowerConeSOCTransform::apply(const Problem& problem) {
  std::vector<Expression> constraints;
  Expression linear_objective = transform_expression(
      problem.objective, &constraints);
  for (const Expression& constr : problem.constraints) {
    constraints.push_back(
        transform_expression(constr, &constraints));
  }
  return {problem.sense, linear_objective, constraints};
}

// Return true if we have a transform for this expression itself
bool have_transform_self(const Expression& expr) {
  if (is_leaf(expr) || is_linear(expr) || is_constraint(expr))
    return true;

  auto iter = kTransforms.find(expr.type());
  if (iter == kTransforms.end()) {
    VLOG(1) << "No transform for " << format_expression(expr);
    return false;
  }
  return true;
}

// Return true if we have a transform this expression and its children
bool have_transform(const Expression& expr) {
  if (!have_transform_self(expr))
    return false;

  for (const Expression& arg : expr.args()) {
    if (!have_transform(arg))
      return false;
  }

  return true;
}

bool PowerConeSOCTransform::accepts(const Problem& problem) {
  for (const Expression& constr : problem.constraints)
    if (!have_transform(constr))
      return false;
  return have_transform(problem.objective);
}
