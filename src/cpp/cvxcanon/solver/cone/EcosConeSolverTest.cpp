
#include "cvxcanon/solver/cone/EcosConeSolver.hpp"
#include "cvxcanon/solver/cone/ConeSolver.hpp"

#include "gtest/gtest.h"

TEST(EcosConeSolverTest, SizeConstraint) {

}

TEST(EcosConeSolverTest, BuildConstraint) {

}

TEST(EcosConeSolverTest, BuildProblem) {

}

TEST(EcosConeSolverTest, SolverStatus) {
ecos_data_->pwork_.info = 0;
EXPECT_EQ(OPTIMAL, get_ecos_status());

ecos_data_->pwork_.info = 1;
EXPECT_EQ(INFEASIBLE, get_ecos_status());

ecos_data_->pwork_.info = 2;
EXPECT_EQ(UNBOUNDED, get_ecos_status());

ecos_data_->pwork_.info = -1;
EXPECT_EQ(USER_LIMT, get_ecos_status());

ecos_data_->pwork_.info = 999;
EXPECT_EQ(ERROR, get_ecos_status());
}

TEST(EcosConeSolverTest, Solve) {

}
