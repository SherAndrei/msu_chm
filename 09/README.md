# Models of Spatial Variables

#### Build
```bash
$ make
```

## Generate points
See help
```bash
$ ./Generate.out
```
To generate 1000 points on a rectange $[0, 2]\times[0, 2]$ using function from `ExactSolution.c` and save the result to `data/input.txt` use
```bash
$ ./Generate.out 1000 2 2 > data/input.txt
```

## Radial Basis Function
See help
```bash
$ ./RBFInterpolator.out --help
```
Using `data/input.txt` calculate values on a rectangle $[0, 2]\times[0, 2]$ with
the grid which consists 200 vertical and 20 horizontal lines and save to `data/output.txt`
```bash
$ ./RBFInterpolator.out 2 2 --Nx 20 --Ny 20 < data/input.txt > data/output.txt
```
To plot solution into `data/solution.png` of a known function $2\sin(10\sqrt{(x-0.5)^2+(y-0.5)^2})$ use gnuplot notation of a function in the next command:
```bash
$ ./print_solution.sh data/input.txt data/output.txt "2*sin(10*sqrt((x-0.5)**2+(y-0.5)**2))" data/solution.png
```

![solution.png](data/solution.png)

If amount of input data is enormous there can be a problem of memory deficiency, i.e.:
```bash
$ ./Generate.out 100000 5 5 > gen.txt
$ ./RBFInterpolator.out 5 5 < gen.txt > out.txt
error: not enough memory for spatial interpolation
```
You can work around this by choosing how many nearest data points participates in the computing of the interpolant at each evaluation by specifiyng `k_neighbors`
```bash
$ ./RBFInterpolator.out 5 5 -k 100 < gen.txt > out.txt
$ ./print_solution.sh gen.txt out.txt "2*sin(10*sqrt((x-0.5)**2+(y-0.5)**2))" data/knn_solution.png
```

![knn_solution.png](data/knn_solution.png)
