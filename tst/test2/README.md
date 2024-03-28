
The reactive Sod Shock tube is a variation of the Sod's Tube case but with a
 mixture of hydrogen-air  with a 2:1:7 mole ratio (CHECK)

The solution is plotted at t=0.2

tested       |      grid     | comment
:----------- |:-------------:| -----------:
gcc 11.4(Linux), 13.x (Mac)       | **128**        |  N-S 1500 steps, no-AMR



After compiling using `make` to run type
```
$ ./main1d.gnu.ex inputs
```
It should run fast (3-5 secs, depending on machine)
The results can be seen by 

```
$ python plt.py
```

![test2plot](images/test2_ref.png)


 note: a slightly different mechanism is used in ref



