#ifndef CVXCANON_TRANSFORM_POWER_CONE_SOC_TRANSFORM_H
#define CVXCANON_TRANSFORM_POWER_CONE_SOC_TRANSFORM_H

#include "cvxcanon/transform/ProblemTransform.hpp"

class PowerConeSOCTransform : public ProblemTransform {
   public:
      PowerConeSOCTransform();
      ~PowerConeSOCTransform();

      bool accepts(const Problem& problem) override;
      Problem apply(const Problem& problem) override;
   private:
      std::vector<Expression> gm_constrs(
         const Expression& t, 
         std::vector<Expression>& expr, 
         std::vector<std::pair<double, double>> p);
};

#endif  // CVXCANON_TRANSFORM_POWER_CONE_SOC_TRANSFORM_H
