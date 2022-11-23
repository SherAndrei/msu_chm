build
```
make release
```
To see the plot of the task
$$
\begin{cases}
y'_1(x)&=y_0(x) \\
y'_0(x)&=-y_1(x)
\end{cases}
$$
with exact solution
$$
\begin{cases}
y_0(x)&=\sin(x) \\
y_1(x)&=\cos(x)
\end{cases}
$$
use
```bash
$ ./sin_and_cos.out 1e-2 out.txt
$ gnuplot
gnuplot> plot 'out.txt' u 1:2 w lp, 'out.txt' u 1:3 w lp, 'out.txt' u 1:5 w lp, 'out.txt' u 1:6 w lp
```
