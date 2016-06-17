#include "LinearExpression.hpp"

#include <unordered_map>
#include <map>
#include <vector>

#include "cvxcanon/expression/ExpressionShape.hpp"
#include "cvxcanon/expression/ExpressionUtil.hpp"
#include "cvxcanon/expression/TextFormat.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include "glog/logging.h"

bool is_constant(const CoeffMap& coeffs) {
  return (coeffs.find(kConstCoefficientId) != coeffs.end() &&
          coeffs.size() == 1);
}

// TODO(mwytock): Replace multiply_matrix_left()/multiply_matrix_right() with a
// kron expressions.
SparseMatrix multiply_matrix_left(const SparseMatrix& A, int n) {
  SparseMatrix coeffs(A.rows()*n, A.cols()*n);
  std::vector<Triplet> tripletList;
  tripletList.reserve(n * A.nonZeros());
  for (int j = 0; j < n; j++) {
    const int row_start = A.rows() * j;
    const int col_start = A.cols() * j;
    for (int i = 0; i < A.outerSize(); i++) {
      for (Matrix::InnerIterator it(A, i); it; ++it) {
        const int row_idx = row_start + it.row();
        const int col_idx = col_start + it.col();
        const double val = it.value();
        tripletList.push_back(Triplet(row_idx, col_idx, val));
      }
    }
  }
  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return coeffs;
}

SparseMatrix multiply_matrix_right(const SparseMatrix& A, int m) {
  SparseMatrix coeffs(m*A.cols(),  m*A.rows());
  std::vector<Triplet> tripletList;
  tripletList.reserve(m * A.nonZeros());
  for (int i = 0; i < A.outerSize(); i++) {
    for (Matrix::InnerIterator it(A, i); it; ++it) {
      const double val = it.value();
      const int row_start = it.col() * m;
      const int col_start = it.row() * m;
      for (int j = 0; j < m; j++) {
        const int row_idx = row_start + j;
        const int col_idx = col_start + j;
        tripletList.push_back(Triplet(row_idx, col_idx, val));
      }
    }
  }
  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return coeffs;
}

std::vector<SparseMatrix> get_add_coefficients(const Expression& expr) {
  std::vector<SparseMatrix> coeffs;
  for (const Expression& arg : expr.args()) {
    // Handle promotion
    coeffs.push_back(
        dim(arg) == 1 ? ones_matrix(dim(expr), 1) : identity(dim(expr)));
  }
  return coeffs;
}

std::vector<SparseMatrix> get_neg_coefficients(const Expression& expr) {
  return {scalar_matrix(-1, dim(expr))};
}

std::vector<SparseMatrix> get_sum_entries_coefficients(const Expression& expr) {
  const int axis = expr.attr<SumEntriesAttributes>().axis;
  const Expression& x = expr.arg(0);

  if (axis == kNoAxis) {
    return {ones_matrix(1, dim(x))};
  } else if (axis == 0) {
    SparseMatrix A = ones_matrix(1, size(x).dims[0]);
    return {multiply_matrix_left(A, size(x).dims[1])};
  } else if (axis == 1) {
    SparseMatrix A = ones_matrix(size(x).dims[1], 1);
    return {multiply_matrix_right(A, size(x).dims[0])};
  } else {
    LOG(FATAL) << "invalid axis " << axis;
  }
}

std::vector<SparseMatrix> get_stack_coefficients(
    const Expression& expr, bool vertical) {
  std::vector<SparseMatrix> coeffs;
  int offset = 0;
  Size expr_size = size(expr);
  for (const Expression& arg : expr.args()) {
    Size arg_size = size(arg);

    /* If VERTICAL, columns that are interleaved. Otherwise, they are
       laid out in order. */
    int column_offset;
    int offset_increment;
    if (vertical) {
      column_offset = expr_size.dims[0];
      offset_increment = arg_size.dims[0];
    } else {
      column_offset = arg_size.dims[0];
      offset_increment = dim(arg);
    }

    std::vector<Triplet> arg_coeffs;
    arg_coeffs.reserve(dim(arg));
    for (int i = 0; i < arg_size.dims[0]; i++) {
      for (int j = 0; j < arg_size.dims[1]; j++) {
        int row_idx = i + (j * column_offset) + offset;
        int col_idx = i + (j * arg_size.dims[0]);
        arg_coeffs.push_back(Triplet(row_idx, col_idx, 1));
      }
    }

    coeffs.push_back(sparse_matrix(dim(expr), dim(arg), arg_coeffs));
    offset += offset_increment;
  }
  return coeffs;
}

std::vector<SparseMatrix> get_hstack_coefficients(const Expression& expr) {
  return get_stack_coefficients(expr, false);
}

std::vector<SparseMatrix> get_vstack_coefficients(const Expression& expr) {
  return get_stack_coefficients(expr, true);
}

std::vector<SparseMatrix> get_reshape_coefficients(const Expression& expr) {
  return {identity(dim(expr))};
}

bool slice_done(const Slice& slice, int i) {
  if (slice.step > 0)
    return i >= slice.stop;
  else
    return i <= slice.stop;
}

std::vector<SparseMatrix> get_index_coefficients(const Expression& expr) {
  const Expression& X = expr.arg(0);
  const int m = size(X).dims[0];
  const int n = size(X).dims[1];
  SparseMatrix coeffs(dim(expr), m*n);

  const Slice& row = expr.attr<IndexAttributes>().keys[0];
  const Slice& col = expr.attr<IndexAttributes>().keys[1];
  int k = 0;
  std::vector<Triplet> tripletList;
  for (int j = col.start; !slice_done(col, j); j += col.step) {
    for (int i = row.start ; !slice_done(row, i); i += row.step) {
      tripletList.push_back(Triplet(k++, j*m+i, 1));
    }
  }

  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return {coeffs};
}

std::vector<SparseMatrix> get_diag_mat_coefficients(const Expression& expr) {
  const int rows = size(expr).dims[0];

  SparseMatrix coeffs(rows, rows * rows);
  std::vector<Triplet> tripletList;
  tripletList.reserve(rows);
  for (int i = 0; i < rows; i++) {
    // index in the extracted vector
    int row_idx = i;
    // index in the original matrix
    int col_idx = i * rows + i;
    tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
  }

  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return {coeffs};
}

std::vector<SparseMatrix> get_diag_vec_coefficients(const Expression& expr) {
  const int rows = size(expr).dims[0];

  SparseMatrix coeffs(rows * rows, rows);
  std::vector<Triplet> tripletList;
  tripletList.reserve(rows);
  for (int i = 0; i < rows; i++) {
    // index in the diagonal matrix
    int row_idx = i * rows + i;
    // index in the original vector
    int col_idx = i;
    tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
  }
  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return {coeffs};
}

std::vector<SparseMatrix> get_transpose_coefficients(const Expression& expr) {
  const int rows = size(expr).dims[0];
  const int cols = size(expr).dims[1];

  SparseMatrix coeffs(rows * cols, rows * cols);
  std::vector<Triplet> tripletList;
  tripletList.reserve(rows * cols);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      int row_idx = rows * j + i;
      int col_idx = i * cols + j;
      tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
    }
  }
  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return {coeffs};
}

std::vector<SparseMatrix> get_kron_coefficients(
    const Expression& expr, const SparseMatrix& constant) {
  int lh_rows = constant.rows();
  int lh_cols = constant.cols();
  int rh_rows = size(expr.arg(1)).dims[0];
  int rh_cols = size(expr.arg(1)).dims[1];

  int rows = rh_rows * rh_cols * lh_rows * lh_cols;
  int cols = rh_rows * rh_cols;
  SparseMatrix coeffs(rows, cols);

  std::vector<Triplet> tripletList;
  tripletList.reserve(rh_rows * rh_cols * constant.nonZeros());
  for ( int k = 0; k < constant.outerSize(); ++k ) {
    for ( Matrix::InnerIterator it(constant, k); it; ++it ) {
      int row = (rh_rows * rh_cols * (lh_rows * it.col())) +
                (it.row() * rh_rows);
      int col = 0;
      for (int j = 0; j < rh_cols; j++) {
        for (int i = 0; i < rh_rows; i++) {
          tripletList.push_back(Triplet(row + i, col, it.value()));
          col++;
        }
        row += lh_rows * rh_rows;
      }
    }
  }

  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return {coeffs};
}

std::vector<SparseMatrix> get_trace_coefficients(const Expression& expr) {
  const int rows = size(expr.arg(0)).dims[0];
  SparseMatrix coeffs(1, rows * rows);
  for (int i = 0; i < rows; i++) {
    coeffs.insert(0, i * rows + i) = 1;
  }
  coeffs.makeCompressed();
  return {coeffs};
}

std::vector<SparseMatrix> get_upper_tri_coefficients(const Expression& expr) {
  const Expression& x = expr.arg(0);
  const int rows = size(x).dims[0];
  const int cols = size(x).dims[1];

  const int entries = size(expr).dims[0];
  SparseMatrix coeffs(entries, rows * cols);

  std::vector<Triplet> tripletList;
  tripletList.reserve(entries);
  int count = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      if (j > i) {
        // index in the extracted vector
        int row_idx = count;
        count++;
        // index in the original matrix
        int col_idx = j * rows + i;
        tripletList.push_back(Triplet(row_idx, col_idx, 1.0));
      }
    }
  }
  coeffs.setFromTriplets(tripletList.begin(), tripletList.end());
  coeffs.makeCompressed();
  return {coeffs};
}

typedef std::vector<SparseMatrix>(*CoefficientFunction)(
    const Expression& expr);

std::unordered_map<int, CoefficientFunction> kCoefficientFunctions = {
  {Expression::ADD, &get_add_coefficients},
  {Expression::DIAG_MAT, &get_diag_mat_coefficients},
  {Expression::DIAG_VEC, &get_diag_vec_coefficients},
  {Expression::HSTACK, &get_hstack_coefficients},
  {Expression::INDEX, &get_index_coefficients},
  {Expression::NEG, &get_neg_coefficients},
  {Expression::RESHAPE, &get_reshape_coefficients},
  {Expression::SUM_ENTRIES, &get_sum_entries_coefficients},
  {Expression::TRACE, &get_trace_coefficients},
  {Expression::TRANSPOSE, &get_transpose_coefficients},
  {Expression::UPPER_TRI, &get_upper_tri_coefficients},
  {Expression::VSTACK, &get_vstack_coefficients},
};

// result += lhs*rhs
void multiply_by_constant(
    const SparseMatrix& lhs, const CoeffMap& rhs_map, CoeffMap* result) {
  for (const auto& iter : rhs_map) {
    const SparseMatrix& rhs = iter.second;
    VLOG(3) << "multiplying\n"
            << "lhs:\n" << matrix_debug_string(lhs)
            << "rhs:\n" << matrix_debug_string(rhs);

    CHECK_EQ(lhs.cols(), rhs.rows());
    SparseMatrix value = lhs*rhs;
    auto ret = result->insert(std::make_pair(iter.first, value));
    if (!ret.second)
      ret.first->second += value;
  }
}

SparseMatrix get_constant(const Expression& expr) {
  CoeffMap coeffs = get_coefficients(expr);
  CHECK(is_constant(coeffs));
  return reshape(coeffs[kConstCoefficientId],
                 size(expr).dims[0], size(expr).dims[1]);
}

CoeffMap get_coefficients(const Expression& expr) {
  VLOG(2) << "get_coefficients\n" << tree_format_expression(expr);

  CoeffMap coeffs;
  if (expr.type() == Expression::CONST || expr.type() == Expression::PARAM) {
    // Special case for constants
    const Constant& constant = (
        expr.type() ==  Expression::CONST ?
        expr.attr<ConstAttributes>().constant :
        expr.attr<ParamAttributes>().constant);
    if (constant.sparse) {
      coeffs[kConstCoefficientId] = reshape(constant.sparse_data, dim(expr), 1);
    } else {
      coeffs[kConstCoefficientId] = to_vector(constant.dense_data).sparseView();
    }
  } else if (expr.type() == Expression::VAR) {
    // Special case for variable
    coeffs[expr.attr<VarAttributes>().id] = identity(dim(expr));
  } else if (expr.type() == Expression::MUL) {
    // Special case for binary mul operator which is guaranteed to have one
    // constant argument by DCP rules.
    CHECK_EQ(expr.args().size(), 2);
    Expression lhs = promote_identity(expr.arg(0), size(expr).dims[0]);
    Expression rhs = promote_identity(expr.arg(1), size(expr).dims[1]);
    CoeffMap lhs_coeffs = get_coefficients(lhs);
    CoeffMap rhs_coeffs = get_coefficients(rhs);

    if (is_constant(lhs_coeffs)) {
      SparseMatrix A = reshape(
          lhs_coeffs[kConstCoefficientId],
          size(lhs).dims[0],
          size(lhs).dims[1]);
      multiply_by_constant(
          multiply_matrix_left(A, size(expr).dims[1]), rhs_coeffs, &coeffs);
    } else if (is_constant(rhs_coeffs)) {
      SparseMatrix A = reshape(
          rhs_coeffs[kConstCoefficientId],
          size(rhs).dims[0],
          size(rhs).dims[1]);
      multiply_by_constant(
          multiply_matrix_right(A, size(expr).dims[0]), lhs_coeffs, &coeffs);
    } else {
      LOG(FATAL) << "multipying two non constants";
    }
  } else if (expr.type() == Expression::KRON) {
    // Special case for binary kron operator, lhs must be constant
    SparseMatrix A = get_constant(expr.arg(0));
    CoeffMap rhs_coeffs = get_coefficients(expr.arg(1));
    std::vector<SparseMatrix> f_coeffs = get_kron_coefficients(expr, A);
    multiply_by_constant(f_coeffs[0], rhs_coeffs, &coeffs);
  } else {
    // General case
    auto iter = kCoefficientFunctions.find(expr.type());
    CHECK(iter != kCoefficientFunctions.end())
        << "no linear coefficients for " << format_expression(expr);
    std::vector<SparseMatrix> f_coeffs = iter->second(expr);
    for (int i = 0; i < expr.args().size(); i++) {
      std::map<int, SparseMatrix> arg_coeffs = get_coefficients(expr.arg(i));
      VLOG(2) << "multiply_by_constant " << format_expression(expr) << " " << i;
      multiply_by_constant(f_coeffs[i], arg_coeffs, &coeffs);
    }
  }

  if (VLOG_IS_ON(2)) {
    VLOG(2) << "get_coefficients done\n" << tree_format_expression(expr);
    for (const auto& iter : coeffs) {
      VLOG(2) <<  "id: " << iter.first << "\n"
              << matrix_debug_string(iter.second);
    }
  }

  return coeffs;
}
