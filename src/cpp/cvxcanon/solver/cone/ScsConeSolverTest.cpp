
#include "cvxcanon/solver/cone/ScsConeSolver.hpp"
#include "cvxcanon/solver/cone/ConeSolver.hpp"

#include <cmath>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

// TEST(ScsConeSolverTest, BuildConstraint) {
//    ConeConstraint con_cone;   

//    int n = 2;
//    int m = 2;

//    std::vector<Triplet> tripletlist;
//    tripletlist.push_back(Triplet(0,0,1));
//    tripletlist.push_back(Triplet(1,0,1));
//    tripletlist.push_back(Triplet(0,1,1));

//    Eigen::SparseMatrix<double, Eigen::RowMajor> A(m,n);
//    A.setFromTriplets(tripletlist.begin(), tripletlist.end());

//    EXPECT_EQ(3, A.nonZeros());

//    DenseVector b(m);
//    EXPECT_EQ(2, b.size());

//    b[0] = 2;
//    b[1] = 1;

//    EXPECT_EQ(2, b[0]);
//    EXPECT_EQ(1, b[1]);

//    std::vector<ConeConstraint> constraints;


//    constraints.push_back({ConeConstraint::FREE, 0, 1});
//    constraints.push_back({ConeConstraint::ZERO, 1, 1});
//    EXPECT_EQ(1,2);
//    int total_size = 0;
//    int sizes[2];

//    std::unordered_map<int, std::vector<ConeConstraint>> constr_map;
//    for (const ConeConstraint constr : constraints) {
//       constr_map[constr.cone].push_back(constr);
//    }

//    EXPECT_EQ(1,2);
//    ScsConeSolver solver;
//    solver.A_coeffs_.clear();
//    solver.b_ = DenseVector(m);
//    solver.num_constrs_ = 0;

//    EXPECT_EQ(1,2);
//    solver.build_scs_constraint(
//                                A, b, constraints, &total_size, sizes);


   // solver.build_scs_constraint(
   //      A, b, constr_map[ConeConstraint::ZERO], &total_size, sizes);

   // EXPECT_EQ(100, total_size);
   
   // solver.build_scs_constraint(
   //      A, b, constr_map[ConeConstraint::NON_NEGATIVE], &total_size, nullptr);

   // EXPECT_EQ(0, total_size);

   // solver.build_scs_constraint(
   //      A, b, constr_map[ConeConstraint::SECOND_ORDER], nullptr, sizes);

   // EXPECT_EQ(0, sizes[0]);

   // solver.build_scs_constraint(
   //      A, b, constr_map[ConeConstraint::SEMIDEFINITE], nullptr, sizes);

   // EXPECT_EQ(0, sizes[0]);

   // solver.build_scs_constraint(
   //      A, b, constr_map[ConeConstraint::PRIMAL_EXPO], &total_size, nullptr);

//    EXPECT_EQ(0, total_size);
// }

TEST(ScsConeSolverTest, Solve) {
   ConeProblem problem;
   ConeConstraint con_cone; 

   int n = 2;
   int m = 2;

   std::vector<Triplet> tripletlist;
   tripletlist.push_back(Triplet(0,0,1));
   tripletlist.push_back(Triplet(1,1,1));

   Eigen::SparseMatrix<double, Eigen::RowMajor> A(m,n);
   A.setFromTriplets(tripletlist.begin(), tripletlist.end());

   problem.A = A;

   EXPECT_EQ(2, problem.A.nonZeros());

   DenseVector b(m);
   EXPECT_EQ(2, b.size());

   b[0] = 2;
   b[1] = 1;

   problem.b = b;

   EXPECT_EQ(2, problem.b[0]);
   EXPECT_EQ(1, problem.b[1]);

   DenseVector c(n);
   EXPECT_EQ(2, c.size());

   c[0] = 1;
   c[1] = 1;

   problem.c = c;

   EXPECT_EQ(1, problem.c[0]);
   EXPECT_EQ(1, problem.c[1]);

   std::vector<ConeConstraint> constraints;

   con_cone.cone = ConeConstraint::ZERO;
   con_cone.offset = 0;
   con_cone.size = 2;

   constraints.push_back({con_cone.cone, con_cone.offset, con_cone.size});

   problem.constraints = constraints;

   ScsConeSolver solver;
   ConeSolution solution = solver.solve(problem);
   EXPECT_EQ(solution.status, OPTIMAL);
   EXPECT_LE(fabs(solution.p_objective_value - 3), .001);
}

TEST(ScsConeSolverTest, SolverStatus) {
   ScsConeSolver solver;

   int status = 1;
   EXPECT_EQ(OPTIMAL, solver.get_scs_status(status));


   status = -2;
   EXPECT_EQ(INFEASIBLE, solver.get_scs_status(status));

   status = -1;
   EXPECT_EQ(UNBOUNDED, solver.get_scs_status(status));

   status = 999;
   EXPECT_EQ(ERROR, solver.get_scs_status(status));
}

