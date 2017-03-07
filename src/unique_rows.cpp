#include "common.h"

// compare the two vectors, return 1 if mat_row > pivot_row, -1 if mat_row < pivot_row and 0 else
int compare_vec(const rowvec& mat_row, const rowvec& pivot_row){
  int v = 0;
  uvec v1 = find(mat_row < pivot_row, 1),
       v2 = find(mat_row > pivot_row, 1);
  if (v1.is_empty() && !v2.is_empty())
    v = -1;
  if (!v1.is_empty() && v2.is_empty())
    v = 1;
  if (!v1.is_empty() && !v2.is_empty())
    v = (v1(0) < v2(0)) ? 1 : -1;
  return v;
}

// quick sort algorithm to sort rows in a matrix
void sortrows(mat& M, const int& left, const int& right){
  if (left < right) {
    int i = left, j = right;
    // find a middle location
    uword mid_loc = (uword) (left+right)/2, pivot_loc = mid_loc;
    // use median of central 5 rows to get a pivot row to split rows
    if (right - left > 5) {
      uvec sortIndex = stable_sort_index(M.col(0).subvec(mid_loc-2, mid_loc+2));
      pivot_loc = as_scalar(find(sortIndex == 2)) + mid_loc - 1;
    }
    // get pivot row
    rowvec pivot_row = M.row(pivot_loc);
    // use pivot row to split matrix into two parts
    while (i <= j) {
      // move i to right if left row > pivot low
      while (compare_vec(M.row( (uword) i), pivot_row) == 1)
        ++i;
      // move j to left if right row < pivot low
      while (compare_vec(M.row( (uword) j), pivot_row) == -1)
        --j;
      // exchange two rows
      if (i <= j) {
        M.swap_rows((uword) i, (uword) j);
        ++i;
        --j;
      }
    }
    // sort right part
    if (j > 0)
      sortrows(M, left, j);
    // sort left part
    if (i < (int) M.n_rows - 1)
      sortrows(M, i, right);
  }
}

//' unique rows
//'
//' Return the unique rows with quick sort algorithm.
//'
//' @param x A matrix.
//' @return A matrix with unique rows.
//' @examples
//' unique_rows(matrix(c(1,1,2,2,3,3,1,1,1,1,4,5),, 2))
//' #      [,1] [,2]
//' # [,1]    1    1
//' # [,2]    2    1
//' # [,3]    3    4
//' # [,4]    3    5
//' @export
// [[Rcpp::export]]
arma::mat unique_rows(arma::mat x){
  chk_mat(x, "x");

  if (x.n_rows > 1) {
    // sort rows
    sortrows(x, 0, x.n_rows - 1);
    // find the unique indecies
    uvec uniIdx = join_cols(ones<uvec>(1), any(x.rows(0, x.n_rows-2) != x.rows(1, x.n_rows-1), 1));
    return x.rows(find(uniIdx));
  } else {
    return x;
  }
}
