
This is test is a 2-D laminar periodic channel flow.
The set-up is based on HAMISH validation case

https://www.ukctrf.com/index.php/benchmarking-of-the-new-software/

with a bulk velocity of 32 m/s

Tested       |          grid | comment
:----------- |:-------------:| -----------:
gcc 14.2     | **8x64**(2 levels)        | 


The flow is forced by  a graidient

F=5e-5 

h=1 (height of channel

$$$
u_{max} = \frac{h^2 F}{8 \mu}
$$


Re=50





 ```python ./plot.py``` 






