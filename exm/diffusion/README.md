
This test check AMR behavipur in a temperature diffusion problem.

The set-up is based on HAMISH

https://www.ukctrf.com/index.php/benchmarking-of-the-new-software/



Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 14.2     | **64x64**(2 levels)        |   200 steps (t=0.02), AMR



### Results
The Initial conditions in temperature looks like

![shock](images/exm_diff_0.png)

At t=0.02 without AMR
![shock](images/exm_diff_noAMR.png)
and with
![shock](images/exm_diff_AMR.png)


The Temeprature distribution at x=0, using ```python ./plot.py``` 


![shock](images/exm_diff_T.png)


