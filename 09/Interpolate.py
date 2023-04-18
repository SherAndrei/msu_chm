#!/usr/bin/env python3

import scipy as sc;
import numpy as np;

import matplotlib.pyplot as plt

def ff(x, y):
  return x**2 + y**2

if __name__ == '__main__':
  x1_list = []
  x2_list = []
  f_list = []
  N = int(input())
  while N > 0:
    x1, x2, f = map(float, input().split())
    x1_list.append(x1)
    x2_list.append(x2)
    f_list.append(f)
    N = N - 1
    
  x1_array = np.array(x1_list)
  x2_array = np.array(x2_list)
  # data = np.array(f_list)
  
  xg, yg = np.meshgrid(x1_array, x2_array, indexing='ij')
  data = ff(xg, yg)
  print(data)  
  
  # interp = sc.RegularGridInterpolator((x, y), data,
  #                               bounds_error=False, fill_value=None)
  
  # interp = sc.RegularGridInterpolator((x1_array, x2_array), data,
  #                                bounds_error=False, fill_value=None)
  # fig = plt.figure()
  # ax = fig.add_subplot(projection='3d')
  # ax.scatter(x1_array.ravel(), x2_array.ravel(), data.ravel(),
  #          s=60, c='k', label='data')
  # print(x1_array)
  # print(x2_array)
  # print(x1g[0])
  # print(x2g[0])
  # print(x1g, x2g)
  # values_array = np.array(f_list)
  

