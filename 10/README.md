### Task (Newton.pdf 1b)

Find a numerical solution of the equation f(x)=0 with an tinitial approximation `x0` and with an accuracy of `eps` using modified Newton's method: replacing the exact calculation of the value derivative to an approximate analog $f'(x)=\frac{f(x+h)-f(x-h)}{2h} + O(h^2)$.

#### Build
```bash
$ make
```

#### Integral
See help
```bash
$ ./ModifiedNewtonRaphson.out
```
To find a numerical solution for $f(x)=0$, $f(x)=x^3 cos^2(x)$ with initial approximation $0.1$, precision $1e-4$ and $h=10^{-8}$:
```bash
$ ./ModifiedNewtonRaphson.out 0.1 1e-4 1e-8
root:   2.269937e-04
```
To show diagnostic of convergence rate with known root $z = 0$
```bash
$ ./ModifiedNewtonRaphson.out 0.1 1e-4 1e-8 0
step    x0              f(x0)           df(x0)          fabs((x0-z)/(x1-z))
1       1.000000e-01    9.900333e-04    2.950233e-02    1.505068e+00
2       6.644220e-02    2.920205e-04    1.314646e-02    1.502220e+00
3       4.422933e-02    8.635375e-05    5.849584e-03    1.500981e+00
4       2.946695e-02    2.556399e-05    2.601136e-03    1.500435e+00
5       1.963895e-02    7.571588e-06    1.156321e-03    1.500193e+00
6       1.309095e-02    2.243049e-06    5.139718e-04    1.500086e+00
7       8.726799e-03    6.645564e-07    2.284421e-04    1.500038e+00
8       5.817718e-03    1.968989e-07    1.015318e-04    1.500017e+00
9       3.878435e-03    5.833955e-08    4.512565e-05    1.500008e+00
10      2.585610e-03    1.728568e-08    2.005592e-05    1.500003e+00
11      1.723736e-03    5.121667e-09    8.913758e-06    1.500001e+00
12      1.149157e-03    1.517529e-09    3.961673e-06    1.500001e+00
13      7.661040e-04    4.496379e-10    1.760744e-06    1.500000e+00
14      5.107359e-04    1.332260e-10    7.825531e-07    1.500000e+00
15      3.404906e-04    3.947437e-11    3.478014e-07    1.500000e+00
root:   2.269937e-04
```
To change $f(x)$ to $\sin(x) - 1$ modify `Functions.c` in the following way
  ```diff
  - double f(double x) { return pow(x, 3.) * pow(cos(x), 2); }
  + double f(double x) { return sin(x) - 1.; }
  ```
