### Task (Newton.pdf 1a)

Find a numerical solution of the equation f(x)=0 with an tinitial approximation `x0` and with an accuracy of `eps` using Newton's method.

#### Build
```bash
$ make
```

#### Integral
See help
```bash
$ ./NewtonRaphson.out
```
To find a numerical solution for $f(x)=0$, $f(x)=(x-1)^2$ with initial approximation $1.0001$ and precision $1e-7$:
```bash
$ ./NewtonRaphson.out 1.0001 1e-7
root:   1.000000e+00
```
To show diagnostic of convergence rate with known root $z = 1$
```bash
$ ./NewtonRaphson.out 1.0001 1e-7 1
step    x0              f(x0)           df(x0)          fabs((x0-z)/(x1-z))
1       1.000100e+00    1.000000e-08    2.000000e-04    2.000000e+00
2       1.000050e+00    2.500000e-09    1.000000e-04    2.000000e+00
3       1.000025e+00    6.250000e-10    5.000000e-05    2.000000e+00
4       1.000012e+00    1.562500e-10    2.500000e-05    2.000000e+00
5       1.000006e+00    3.906250e-11    1.250000e-05    2.000000e+00
6       1.000003e+00    9.765625e-12    6.250000e-06    2.000000e+00
7       1.000002e+00    2.441406e-12    3.125000e-06    2.000000e+00
8       1.000001e+00    6.103516e-13    1.562500e-06    2.000000e+00
9       1.000000e+00    1.525879e-13    7.812500e-07    2.000000e+00
root:   1.000000e+00
```
To change $f(x)$ to $\sin(x) - 1$ modify `Functions.c` in the following way
  ```diff
  - double f(double x) { return pow(x - 1., 2); }
  + double f(double x) { return sin(x) - 1.; }
  - double df(double x) { return 2. * (x - 1.); }
  + double df(double x) { return cos(x); }
  ```
