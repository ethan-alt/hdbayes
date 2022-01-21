
// find index of x, j, such that x0 is closest to x[j] without
// going over. Uses binary search algorithm
int findClosestIndex(real x0, vector x) {
  int K = rows(x);
  int i = 1;
  int j = K;
  int mid;
  // check corner cases
  if ( x0 < x[2] )
    return 1;
  if ( x0 == x[K] ) {
    return K;
  }
  // conduct binary search
  while ( i <= j ) {
    mid = (i + j) / 2;
    // if x0 < x[mid], index must lie in left half
    if ( x0 < x[mid] ) {
      // if x0 is larger than x[mid-1], return mid-1; else update j
      if ( mid > 2 &&  x0 > x[mid - 1] )
        return mid - 1;
      j = mid;
    }
    // otherwise, index must lie in right half
    else {
      // if x0 is less than x[mid + 1], return mid; else update i
      if ( mid < K && x0 < x[mid + 1] )
        return mid;
      i = mid + 1;
    }
  }
  reject("Error in finding midpoint");
  return(0); // never reached
}

// approximate lognc of power prior
//
// * @param a0       power prior param to obtain lognc
// * @param a0vec    fine grid of power prior parameters for which we have estimates
// * @param lognca0  estimate lognc pertaining to fine grid a0vec
//
// * @return linearly interpolated log normalizing constant.
real pp_lognc(real a0, vector a0vec, vector lognca0) {
  // find index of a0vec closest to a0
  int i = findClosestIndex(a0, a0vec);
  // if not exact match, use linear interpolation to get estimated lognc
  if ( a0 != a0vec[i] ) {
    real x1 = a0vec[i];
    real x2 = a0vec[i + 1];
    real y1 = lognca0[i];
    real y2 = lognca0[i + 1];
    return y1 + (y2 - y1) * (a0 - x1) / (x2 - x1);
  }
  return lognca0[i];
}
