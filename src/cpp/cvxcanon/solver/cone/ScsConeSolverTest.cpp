
#include "cvxcanon/solver/cone/ScsConeSolver.hpp"
#include "cvxcanon/solver/cone/ConeSolver.hpp"

#include "gtest/gtest.h"

ConeProblem problem;

//TODO(fabioftv): Input problem.A, problem.b, problem.c, and ConeConstraints 

TEST(ScsConeSolverTest, BuildConstraint) {

}

TEST(ScsConeSolverTest, BuildProblem) {

}

TEST(ScsConeSolverTest, SolverStatus) {
scs_data_->info_.status = "Solved";
EXPECT_EQ(OPTIMAL, get_scs_status());

scs_data_->info_.status = "Infeasible";
EXPECT_EQ(INFEASIBLE, get_scs_status());

scs_data_->info_.status = "Unbounded";
EXPECT_EQ(UNBOUNDED, get_scs_status());

scs_data_->info_.status = "Test";
EXPECT_EQ(ERROR, get_scs_status());
}

TEST(ScsConeSolverTest, Solve) {

}
