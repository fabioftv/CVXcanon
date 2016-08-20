#ifndef CVXCANON_TRANSFORM_POWER_CONE_TRANSFORM_H
#define CVXCANON_TRANSFORM_POWER_CONE_TRANSFORM_H

#include "cvxcanon/transform/ProblemTransform.hpp"

class PowerConeTransform : public ProblemTransform {
 public:
  bool accepts(const Problem& problem) override;
  Problem apply(const Problem& problem) override;
};

#endif  // CVXCANON_TRANSFORM_POWER_CONE_TRANSFORM_H
