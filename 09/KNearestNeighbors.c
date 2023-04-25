#include "KNearestNeighbors.h"

#include <stdlib.h>

static inline int Less(const DistanceAndIndex *p_lhs, const DistanceAndIndex *p_rhs) {
  const double lhs = p_lhs->distance;
  const double rhs = p_rhs->distance;
  return lhs < rhs;
}

static inline void IterSwap(DistanceAndIndex *lhs, DistanceAndIndex *rhs) {
  DistanceAndIndex tmp = *lhs;
  *lhs = *rhs;
  *rhs = tmp;
}

static inline int Partition(DistanceAndIndex *values, int start, int pivot, int end) {
  int i = start, j = end;

  while (i < j) {
    while (Less(values + i, values + pivot))
      i++;
    while (Less(values + pivot, values + j))
      j--;

    if (i < j) {
      IterSwap(values + i, values + j);
      if (!Less(values + i, values + j))
        --j;
    }
  }
  return i;
}

static inline void QuickSelect(DistanceAndIndex *values, int k, int start, int end) {
  int pivot = Partition(values, start, (start + end) / 2, end);
  int length = pivot - start + 1;
  if (k == length) return;
  if (k < length)
    QuickSelect(values, k, start, pivot - 1);
  else
    QuickSelect(values, k - length, pivot + 1, end);
}

void KNearestNeighbors(const Point *all_points, int N, Point target, int k,
                       DistanceAndIndex *distances) {
  for (int i = 0; i < N; ++i) {
    distances[i].distance = Distance(target, all_points[i]);
    distances[i].index = i;
  }

  QuickSelect(distances, k, 0, N - 1);
}
