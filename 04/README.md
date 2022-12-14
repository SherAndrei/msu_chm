Build with
```bash
$ make release
```
To see the result of the Fourier method for N use
```bash
$ ./Fouirer.out N out.txt
```
To compare the result of the method and the exact solution use
```bash
$ gnuplot -c compare.gnuplot out.txt cmp.png
```
or
To get text file with different N and corresponding errors use
```bash
$ ./collect_errors.sh > err.txt
```
To plot logariphmic scale of errors into .png file use
```bash
$ gnuplot -c scale.gnuplot err.txt scale.png
```
