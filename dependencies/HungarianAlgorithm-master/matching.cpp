#include "Utils.h"
#include "matching.h"

/*
 * reduce
 * reduces matrix based on row and column minimums
 */
void reduce(MatrixReal& m) {
  // subtract row minimum from each row
  for (int i=0; i<m.rows(); i++) {
    Real minElement = m.row(i).minCoeff();
    m.row(i).array() -= minElement;
  }
}

/*
 * hasMark
 * if there is a starred/primed zero in the given row/col, returns it's index
 * else, returns -1
 */
inline int hasMark(MatrixReal::RowXpr& v) {
    for (int i=0; i<v.size(); i++) {
        if (v(i)) {
            return i;
        }
    }
    return -1;
}
inline int hasMark(MatrixReal::ColXpr& v) {
    for (int i = 0; i < v.size(); i++) {
        if (v(i)) {
            return i;
        }
    }
    return -1;
}

/*
 * swapStarsAndPrimes
 * Swap stars and primes based on step 5 of Hungarian algorithm
 * Z0 is uncovered primed zero we've found
 * Z1 is the stared zero in the column of Z0 (if any)
 * Z2 is the primed zero in the row of Z1 (will always be one)
 * ...continue series until we reach a primed zero with no starred zero in its column
 * Unstar each starred zero, star each primed zero, erase all primes and uncover every line in the matrix
 */
void swapStarsAndPrimes(int i, int j, MatrixReal& stars, MatrixReal& primes) {
  int primeRow = i;
  int primeCol = j;
  
  bool done = false;
  while (!done) {
    // find row index of row that has a 0* in the same col as the current 0'
    int starInPrimeColRow = hasMark(stars.col(primeCol));
    
    if (starInPrimeColRow < 0) {
      // star the prime we're looking at
      primes(primeRow, primeCol) = 0;
      stars(primeRow, primeCol) = 1;
      done = true;
    }
    else {
      // find which col has a 0' in the same row as z1
      int primeInStarRowCol = hasMark(primes.row(starInPrimeColRow));
      
      // star first primed zero
      primes(primeRow, primeCol) = 0;
      stars(primeRow, primeCol) = 1;
      //primes(starInPrimeColRow, primeInStarRowCol) = 0;
      //stars(starInPrimeColRow, primeInStarRowCol) = 1;
      
      // unstar starred zero
      stars(starInPrimeColRow, primeCol) = 0;
      
      // set index of last prime, will check it's column for 0*s next
      primeRow = starInPrimeColRow;
      primeCol = primeInStarRowCol;
    }
  }
  // clear primes
  primes.fill(0);
}

/*
 * findMatching
 * implementation of the Hungarian matching algorithm
 * referenced from: http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
 */
bool HungarianMatching(const MatrixReal& m, MatrixReal& result, EHungarianMatchingType type, int maxsteps) {
  MatrixReal n = m; // make a copy of m for reducing
  int dim = n.rows(); // dimension of matrix, used for checking if we've reduced
                      // the matrix enough yet
  
  MatrixReal stars(m.rows(), m.cols()); // matrix for storing our "starred" 0s (0*)
  stars.setZero();
  MatrixReal primes(m.rows(), m.cols()); // matrix for storing our "primed" 0s (0')
  primes.setZero();
  VectorReal rowCover(m.rows()); // keep track of which rows are "covered"
  rowCover.setZero();
  VectorReal colCover(m.cols()); // keep track of which columns are "covered"
  colCover.setZero();
  
  // to do maximization rather than minimization, we have to
  // transform the matrix by subtracting every value from the maximum
  if (type == HUNGARIAN_MATCH_MAX) {
    Real max = n.maxCoeff();
    MatrixReal maxMat(n.rows(), n.cols());
    maxMat.fill(max);
    n = maxMat - n;
  }
  
  // Step 1 
  // Reduce matrix
  reduce(n);
  
  // Step 2
  // Find a zero in the matrix. If there is no starred zero in 
  // its row or column, star Z. Repeat for each element in the matrix.
  for (int i=0; i<n.rows(); i++) {
    for (int j=0; j<n.cols(); j++) {
      if (n(i,j) == 0 && !rowCover(i) && !colCover(j)) {
        stars(i,j) = 1;
        rowCover(i) = 1;
        colCover(j) = 1;
      }
    }
  }
  // covers need to be cleared for following steps
  rowCover.fill(0);
  colCover.fill(0);
  
  int step = 0;
  while (true) {
    // Step 3
    // Cover all columns that have a starred zero
    // If the number of columns with starred zeroes equals the matrix
    // dimensions, we are done! Otherwise, move on to step 4.
    step3:
    for (int j=0; j<n.cols(); j++) {
      if (hasMark(stars.col(j)) >= 0) {
        colCover(j) = 1;
      }
    }
    if (colCover.sum() == dim) {
      result = stars;
      return true;
    }
    
    // Step 4
    // Find a non-covered zero and prime it
    step4:
    for (int i=0; i<n.rows(); i++) {
        if (rowCover(i)) {
            continue;
        }
        for (int j=0; j<n.cols(); j++) {
            if (colCover(j)) {
                continue;
            }
            if (n(i,j) == 0) {
                primes(i,j) = 1;
                // if no starred zero in the row...
                MatrixReal::RowXpr row = stars.row(i);
                int mark = hasMark(row);
                if (mark < 0) {
                    // Step 5
                    // swap stars and primes            
                    swapStarsAndPrimes(i, j, stars, primes);
    
                    // clear lines
                    rowCover.setZero();
                    colCover.setZero();
            
                    goto step3;
                } else {
                    // cover row
                    rowCover(i) = 1;
            
                    // uncover column of the starred zero in the same row
                    colCover(mark) = 0;
                }
            }
        }
    }
    
    // Step 6
    // Should now be no more uncovered zeroes
    // Get the minimum uncovered element
	Real min = std::numeric_limits<Real>::infinity();
    for (int i=0; i<n.rows(); i++) {
        if (rowCover(i)) {
            continue;
        }
        for (int j=0; j<n.cols(); j++) {
            if (colCover(j)) {
                continue;
            }
            if (n(i,j) < min) {
                min = n(i,j);
            }
        }
    }
    
    // Subtract minimum from uncovered elements, add it to elements covered twice
    for (int i=0; i<n.rows(); i++) {
        if (rowCover(i)) {
            for (int j = 0; j < n.cols(); j++) {
                if (colCover(j)) {
                    n(i, j) += min;
                }
            }
        } else {
            for (int j = 0; j < n.cols(); j++) {
                if (!colCover(j)) {
                    n(i, j) -= min;
                }
            }
        }
    }
    
	if (++step == maxsteps) {
		result = n;
		return false;
	}
    goto step4;
  }

  return true;
}
