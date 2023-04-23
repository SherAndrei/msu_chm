#include "KNearestNeighbors.h"

#include <stdlib.h>

static inline int ByDistance(const DistanceAndIndex *p_lhs, const DistanceAndIndex *p_rhs) {
  const double lhs = p_lhs->distance;
  const double rhs = p_rhs->distance;
  return (lhs > rhs) - (lhs < rhs);
}

static inline void IterSwap(DistanceAndIndex *lhs, DistanceAndIndex *rhs) {
  DistanceAndIndex tmp = *lhs;
  *lhs = *rhs;
  *rhs = tmp;
}

static inline int Pivot(DistanceAndIndex *values, int pInd, int end) {
  int i = 0, j = end;

  while (i <= j) {
    while (ByDistance(values + i, values + pInd) < 0 && i <= j) {
      i++;
    }
    while (ByDistance(values + j, values + pInd) > 0 && j >= i) {
      j--;
    }

    if (i == j)
      break;
    else if (ByDistance(values + i, values + j) == 0) {
      i++;
      continue;
    }

    IterSwap(values + i, values + j);
  }
  return j;
}

static inline void PartialQuicksort(DistanceAndIndex *values, int k, int start, int end) {
  int piv = Pivot(values, (start + end) / 2, end);
  int length = piv - start + 1;
  if (k < length)
    PartialQuicksort(values, k, start, piv - 1);
  else if (k > length)
    PartialQuicksort(values, k - length, piv + 1, end);
  else
    return;
}

void KNearestNeighbors(const Point *all_points, int N, Point target, int k,
                       DistanceAndIndex *distances) {
  for (int i = 0; i < N; ++i) {
    distances[i].distance = Distance(target, all_points[i]);
    distances[i].index = i;
  }

  PartialQuicksort(distances, k, 0, N - 1);
}
