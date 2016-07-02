
#include "cvxcanon/solver/cone/EcosConeSolver.hpp"
#include "cvxcanon/solver/cone/ConeSolver.hpp"

#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

TEST(EcosConeSolverTest, SizeConstraint) {
   std::vector<ConeConstraint> constraints;
/*
   constraints.push_back(ConeConstraint::ZERO);
   constraints.push_back(ConeConstraint::ZERO);
   constraints.push_back(ConeConstraint::ZERO);
   constraints.push_back(ConeConstraint::SECOND_ORDER);
   constraints.push_back(ConeConstraint::NON_NEGATIVE);
   constraints.push_back(ConeConstraint::ZERO);
   constraints.push_back(ConeConstraint::PRIMAL_EXPO);
   constraints.push_back(ConeConstraint::SECOND_ORDER);
   constraints.push_back(ConeConstraint::ZERO);
   constraints.push_back(ConeConstraint::NON_NEGATIVE);
*/

   std::unordered_map<int, std::vector<ConeConstraint>> constr_map;
   for (const ConeConstraint constr : constraints) {
      constr_map[constr.cone].push_back(constr);
   }

   int size = 0;
   EcosConeSolver solver;
   solver.define_size_ecos_constraint(constr_map[ConeConstraint::ZERO], size);
   EXPECT_EQ(0, size);

   solver.define_size_ecos_constraint(constr_map[ConeConstraint::NON_NEGATIVE],
      size);
   EXPECT_EQ(0, size);

   solver.define_size_ecos_constraint(constr_map[ConeConstraint::SECOND_ORDER],
      size);
   EXPECT_EQ(0, size);

   solver.define_size_ecos_constraint(constr_map[ConeConstraint::PRIMAL_EXPO],
      size);
   EXPECT_EQ(0, size);
}

TEST(EcosConeSolverTest, BuildConstraint) {
   int n = 2;
   int m = 4;

   std::vector<Triplet> tripletlistA;
   tripletlistA.push_back(Triplet(0,0,1));
   tripletlistA.push_back(Triplet(1,0,1));
   tripletlistA.push_back(Triplet(3,0,1));
   tripletlistA.push_back(Triplet(0,1,1));
   tripletlistA.push_back(Triplet(2,1,1));
   tripletlistA.push_back(Triplet(3,1,-1));

   Eigen::SparseMatrix<double, Eigen::RowMajor> A(m,n);
   A.setFromTriplets(tripletlistA.begin(), tripletlistA.end());

   EXPECT_EQ(6, A.nonZeros());

   DenseVector b(m);
   EXPECT_EQ(4, b.size());

   b[0] = 2;
   b[1] = 1;
   b[2] = 1;
   b[3] = 0;

   EXPECT_EQ(2, b[0]);
   EXPECT_EQ(1, b[1]);
   EXPECT_EQ(1, b[2]);
   EXPECT_EQ(0, b[3]);

   std::vector<ConeConstraint> constraints;

//(fabioftv): Select both constraints below:

//   constraints.push_back(ConeConstraint::ZERO);
//   constraints.push_back(ConeConstraint::ZERO);

//(fabioftv): Select below at most two constraints:

//   constraints.push_back(ConeConstraint::NON_NEGATIVE);
//   constraints.push_back(ConeConstraint::NON_NEGATIVE);
//   constraints.push_back(ConeConstraint::SECOND_ORDER);
//   constraints.push_back(ConeConstraint::SECOND_ORDER);
//   constraints.push_back(ConeConstraint::SEMIDEFINITE);
//   constraints.push_back(ConeConstraint::SEMIDEFINITE);
//   constraints.push_back(ConeConstraint::PRIMAL_EXPO);
//   constraints.push_back(ConeConstraint::PRIMAL_EXPO);

   std::unordered_map<int, std::vector<ConeConstraint>> constr_map;
   for (const ConeConstraint constr : constraints) {
      constr_map[constr.cone].push_back(constr);
   }

   EcosConeSolver solver;
   solver.build_ecos_constraint(A, b, constr_map[ConeConstraint::ZERO]);
   solver.build_ecos_constraint(A, b, constr_map[ConeConstraint::NON_NEGATIVE]);
   solver.build_ecos_constraint(A, b, constr_map[ConeConstraint::SECOND_ORDER]);
   solver.build_ecos_constraint(A, b, constr_map[ConeConstraint::SEMIDEFINITE]);
   solver.build_ecos_constraint(A, b, constr_map[ConeConstraint::PRIMAL_EXPO]);
/*   
   EXPECT_EQ(3, solver.A_coeffs_.size());
   EXPECT_EQ(3, solver.G_coeffs_.size());
*/
}

TEST(EcosConeSolverTest, BuildProblem) {
   ConeProblem problem;
   ConeSolution *solution;

   int n = 2;
   int m = 4;

   std::vector<Triplet> tripletlistA;
   tripletlistA.push_back(Triplet(0,0,1));
   tripletlistA.push_back(Triplet(1,0,1));
   tripletlistA.push_back(Triplet(3,0,1));
   tripletlistA.push_back(Triplet(0,1,1));
   tripletlistA.push_back(Triplet(2,1,1));
   tripletlistA.push_back(Triplet(3,1,-1));

   Eigen::SparseMatrix<double, Eigen::RowMajor> A(m,n);
   A.setFromTriplets(tripletlistA.begin(), tripletlistA.end());

   problem.A = A;   

   EXPECT_EQ(6, problem.A.nonZeros());

   DenseVector b(m);
   EXPECT_EQ(4, b.size());

   b[0] = 2;
   b[1] = 1;
   b[2] = 1;
   b[3] = 0;

   problem.b = b;

   EXPECT_EQ(2, problem.b[0]);
   EXPECT_EQ(1, problem.b[1]);
   EXPECT_EQ(1, problem.b[2]);
   EXPECT_EQ(0, problem.b[3]);

   DenseVector c(n);
   EXPECT_EQ(2, c.size());

   c[0] = 1;
   c[1] = 1;

   problem.c = c;

   EXPECT_EQ(1, problem.c[0]);
   EXPECT_EQ(1, problem.c[1]);

   std::vector<ConeConstraint> constraints;

//(fabioftv): Select both constraints below:

//   constraints.push_back(ConeConstraint::ZERO);
//   constraints.push_back(ConeConstraint::ZERO);

//(fabioftv): Select below at most two constraints:

//   constraints.push_back(ConeConstraint::NON_NEGATIVE);
//   constraints.push_back(ConeConstraint::NON_NEGATIVE);
//   constraints.push_back(ConeConstraint::SECOND_ORDER);
//   constraints.push_back(ConeConstraint::SECOND_ORDER);
//   constraints.push_back(ConeConstraint::SEMIDEFINITE);
//   constraints.push_back(ConeConstraint::SEMIDEFINITE);
//   constraints.push_back(ConeConstraint::PRIMAL_EXPO);
//   constraints.push_back(ConeConstraint::PRIMAL_EXPO);

   problem.constraints = constraints;

   EcosConeSolver solver;
//   solver.build_ecos_problem(problem, solution);
}

TEST(EcosConeSolverTest, SolverStatus) {
EcosConeSolver solver;
int exitflag = 0;
EXPECT_EQ(OPTIMAL, solver.get_ecos_status(exitflag));

exitflag = 1;
EXPECT_EQ(INFEASIBLE, solver.get_ecos_status(exitflag));

exitflag = 2;
EXPECT_EQ(UNBOUNDED, solver.get_ecos_status(exitflag));

exitflag = -1;
EXPECT_EQ(USER_LIMIT, solver.get_ecos_status(exitflag));

exitflag = 999;
EXPECT_EQ(ERROR, solver.get_ecos_status(exitflag));
}

TEST(EcosConeSolverTest, Solve) {
   ConeProblem problem;
   ConeSolution *solution;

   int n = 2;
   int m = 4;

   std::vector<Triplet> tripletlistA;
   tripletlistA.push_back(Triplet(0,0,1));
   tripletlistA.push_back(Triplet(1,0,1));
   tripletlistA.push_back(Triplet(3,0,1));
   tripletlistA.push_back(Triplet(0,1,1));
   tripletlistA.push_back(Triplet(2,1,1));
   tripletlistA.push_back(Triplet(3,1,-1));

   Eigen::SparseMatrix<double, Eigen::RowMajor> A(m,n);
   A.setFromTriplets(tripletlistA.begin(), tripletlistA.end());

   problem.A = A;   

   EXPECT_EQ(6, problem.A.nonZeros());

   DenseVector b(m);
   EXPECT_EQ(4, b.size());

   b[0] = 2;
   b[1] = 1;
   b[2] = 1;
   b[3] = 0;

   problem.b = b;

   EXPECT_EQ(2, problem.b[0]);
   EXPECT_EQ(1, problem.b[1]);
   EXPECT_EQ(1, problem.b[2]);
   EXPECT_EQ(0, problem.b[3]);

   DenseVector c(n);
   EXPECT_EQ(2, c.size());

   c[0] = 1;
   c[1] = 1;

   problem.c = c;

   EXPECT_EQ(1, problem.c[0]);
   EXPECT_EQ(1, problem.c[1]);

   std::vector<ConeConstraint> constraints;

//(fabioftv): Select both constraints below:

//   constraints.push_back(ConeConstraint::ZERO);
//   constraints.push_back(ConeConstraint::ZERO);

//(fabioftv): Select below at most two constraints:

//   constraints.push_back(ConeConstraint::NON_NEGATIVE);
//   constraints.push_back(ConeConstraint::NON_NEGATIVE);
//   constraints.push_back(ConeConstraint::SECOND_ORDER);
//   constraints.push_back(ConeConstraint::SECOND_ORDER);
//   constraints.push_back(ConeConstraint::SEMIDEFINITE);
//   constraints.push_back(ConeConstraint::SEMIDEFINITE);
//   constraints.push_back(ConeConstraint::PRIMAL_EXPO);
//   constraints.push_back(ConeConstraint::PRIMAL_EXPO);

   problem.constraints = constraints;

   EcosConeSolver solver;
//   solver.solve(problem);

//   ConeSolution solution;
//   EXPECT_EQ(2, solution->p_objective_value);
}
