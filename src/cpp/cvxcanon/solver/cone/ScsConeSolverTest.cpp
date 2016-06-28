
#include <string>

#include "gtest/gtest.h"

#include "cvxcanon/solver/cone/ScsConeSolver.hpp"
#include "cvxcanon/solver/cone/ConeSolver.hpp"

// SCS Environment
namespace scs {
#include "scs/include/scs.h"
#include "scs/include/util.h"
#include "scs/linsys/amatrix.h"
typedef double scs_float;
typedef int scs_int;
}  // namespace scs

struct ScsConeSolver::ScsData {
  // SCS data structures
  scs::Data data_;
  scs::Cone cone_;
  scs::Info info_;
  scs::Sol sol_;
  scs::Settings settings_;

  // SCS supporting data structures
  scs::AMatrix A_matrix_;
};

ScsConeSolver::ScsConeSolver() : scs_data_(new ScsData()) {}
ScsConeSolver::~ScsConeSolver() {}

TEST(ScsConeSolverTest, SolverStatus) {

scs_data_->info_.status = "Solved";
EXPECT_EQ("OPTIMAL", SolverStatus());
/*
SolverStatus ScsConeSolver::get_scs_status() {
  if (strcmp(scs_data_->info_.status, "Solved") == 0) {
    return OPTIMAL;
  } else if (strcmp(scs_data_->info_.status, "Infeasible") == 0) {
    return INFEASIBLE;
  } else if (strcmp(scs_data_->info_.status, "Unbounded") == 0) {
    return UNBOUNDED;
  } else {
    return ERROR;
  }
}
*/


}
