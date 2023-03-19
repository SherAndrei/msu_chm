### Task

* For a given rectangle with sides $Lx$, $Ly$ construct its triangulation
* Calculate the approximate value of the 2-dimensional definite integral using next quadrature
$$I(f)=\underset{T}{\int\int}f(x)dx\approx S(f)=\frac{1}{3}S(T)(f(A)+f(B)+f(C))$$
where $T$ is a triangle in the plane, $S(T)$ its area, $A$, $B$, $C$ the midpoints of the sides.

#### Build
```bash
$ make
```

#### Integral
See help
```bash
$ ./Triangulate.out
```
```bash
$ ./Integrate.out --help
```
To construct triangualtion for a square 1 by 1 and save it to `output.txt`:
```
$ ./Triangulate.out 1 1 > output.txt
```
To find out result of the
$$\underset{[0,1] \times [0,1]}{\int\int}(x^4_1 + x^2_1 x^2_2 + x^4_2) dx_1 dx_2=\int_0^1 \left(\frac{x_1^5}{5} + \frac{x_2^2 x_1^3}{3} + x_2^4 x_1\right)_0^1 dx_2 = \int_0^1 \left(\frac{1}{5} + \frac{x_2^2}{3} + x_2^4\right) dx_2 = \frac{23}{45}\approx 0.5111$$
using generated `output.txt`:
```
$ ./Integrate.out < output.txt
```
Using single line
```bash
$ ./Triangulate.out 1 1 | ./Integrate.out
0.51111111120834
```
