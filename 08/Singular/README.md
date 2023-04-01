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
step      x0                      f(x0)                   df(x0)                  fabs(x1-z)
1       1.00010000000000        0.00000001000000        0.00020000000000        0.00005000000000
2       1.00005000000000        0.00000000250000        0.00010000000000        0.00002500000000
3       1.00002500000000        0.00000000062500        0.00005000000000        0.00001250000000
4       1.00001250000000        0.00000000015625        0.00002500000000        0.00000625000000
5       1.00000625000000        0.00000000003906        0.00001250000000        0.00000312500000
6       1.00000312500000        0.00000000000977        0.00000625000000        0.00000156250000
7       1.00000156250000        0.00000000000244        0.00000312500000        0.00000078125000
8       1.00000078125000        0.00000000000061        0.00000156250000        0.00000039062500
9       1.00000039062500        0.00000000000015        0.00000078125000        0.00000019531250
root:   1.000000e+00
```
To change $f(x)$ to $\sin(x) - 1$ modify `Functions.c` in the following way
  ```diff
  - double f(double x) { return pow(x - 1., 2); }
  + double f(double x) { return sin(x) - 1.; }
  - double df(double x) { return 2. * (x - 1.); }
  + double df(double x) { return cos(x); }
  ```
