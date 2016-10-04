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

typedef Expression(*TransformFunction)(
   const Expression& expr,
   std::vector<Expression>* constraints);

Expression transform_geo_mean_ineq(
   const Expression& expr,
   std::vector<Expression>* constraints) {
  
   const double p = expr.attr<GeoMeanIneqAttributes>().p;
   const Expression& x = expr.arg(0);
   const Expression& y = expr.arg(0);
   Expression t = epi_var(expr, "geo_mean_ineq");
   Expression vec_t = reshape(t, dim(t), 1);
   Expression vec_x = reshape(x, dim(x), 1);
   Expression vec_y = reshape(y, dim(y), 1);
   
   constraints->push_back(soc(add(vec_x, vec_y), vstack({add(vec_x, vec_y), 
      mul(constant(2), vec_t)})));

   return t;
}

std::vector<Expression> PowerConeSOCTransform::gm_constrs(const Expression& t, 
   std::vector<Expression>& expr, std::vector<std::pair<double, double>> p) {

   int i, j;
   bool aux = true;
   GeoMeanIneq geo_mean;

   assert (geo_mean.weight_vector_test(p) == true);

   std::vector<std::pair<double, double>> w(p.size());
   w = geo_mean.dyad_completion(p);

   std::vector<std::pair<std::vector<std::pair<double, double>>, 
      std::pair<std::vector<std::pair<double, double>>, 
         std::vector<std::pair<double, double>>>>> tree;
   tree = geo_mean.decompose(w);

//TODO(fabioftv): Check this with Steven.
   std::vector<std::pair<std::vector<std::pair<double, double>>, 
      Expression>> d;
/*
   std::vector<Expression> d;
   for (i = 0; i < expr.size(); i++)
      d.push_back(var(t.args().size(), t.args().size(),
                  t.attr<VarAttributes>().id));
   //d[w] = t ????
*/

   if (expr.size() < w.size()){
      expr.push_back(t);
   }
   
   assert (expr.size() == w.size());

   for (i = 0; i < expr.size(); i++){
      aux = true;
      for(j = 0; j < w.size(); j++){
         if (w[j].first / w[j].second < 0){
            aux = false;
            j = w.size();
         }
      }
      if (aux = true){
         std::vector<std::pair<double, double>> 
            tmp(w.size(), std::make_pair(0, 0));
         tmp[i].first = tmp[i].second = 1;
         d[i].first = tmp;
         d[i].second = expr[i];
      }
   }

   std::vector<Expression> constraints;
//   transform_geo_mean_ineq(expr, constraints);

   return constraints;
}

std::unordered_map<int, TransformFunction> kTransforms = {
  {Expression::GEO_MEAN_INEQ, &transform_geo_mean_ineq},
};

Expression transform_expression(
   const Expression& expr, 
   std::vector<Expression>* constraints) {

   //First, transform the children.
   std::vector<Expression> linear_args;
   for (const Expression& arg : expr.args())
      linear_args.push_back(transform_expression(arg, constraints));

   //Second, a new expresion is formed with linear args.
   Expression output = Expression(expr.type(), linear_args, expr.attr_ptr());

   //Finally, transform non-linear functions, if necessary.
   VLOG(2) << "transform_func: " << format_expression(expr);
   auto iter = kTransforms.find(expr.type());
   if (iter != kTransforms.end())
      output = iter->second(output, constraints);

   return output;
}

Problem PowerConeSOCTransform::apply(const Problem& problem) {
  std::vector<Expression> constraints;
  Expression linear_objective = transform_expression(problem.objective, 
                                                     &constraints);
  for (const Expression& constr : problem.constraints) {
    constraints.push_back(transform_expression(constr, &constraints));
  }
  return {problem.sense, linear_objective, constraints};
}

bool PowerConeSOCTransform::have_transform_self(const Expression& expr) {
   if (is_leaf(expr) || is_linear(expr) || is_constraint(expr))
      return true;
   auto iter = kTransforms.find(expr.type());
   if (iter == kTransforms.end()) {
      VLOG(1) << "No transform for " << format_expression(expr);
      return false;
   }
   return true;
}

bool PowerConeSOCTransform::have_transform(const Expression& expr) {
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

