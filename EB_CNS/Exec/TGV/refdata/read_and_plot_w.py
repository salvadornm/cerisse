#! /usr/bin/env python
# -*- coding: UTF-8 -*-
# Small python script to extract and plot
# the data contained in the file wn_slice_x0_08000.out
# Other example of script:
#  f = fopen(file, "rb") ;
#  for i
#   for j
#     fread(W(i,j), 1, sizeof(double), f)
#   end
#  end


from numpy import *

import pylab as p
from struct import *

Pi = arccos(-1)
x = arange(-Pi,Pi-2.*Pi/512.,2.*Pi/512.)
y = arange(-Pi,Pi-2.*Pi/512.,2.*Pi/512.)

size = (len(x),len(y))

w = zeros(size)

file = open('wn_slice_x0_08000.out', 'rb')
data = file.read()

n = 0
for i in range(len(x)):
  for j in range(len(y)):
    a = unpack('d', data[n:n+8])
    w[j][i] = a[0]
    n = n+8

fig=p.figure()
p.contour(x,y,w, levels = [1,5,10,20,30])
p.axis('normal')
p.colorbar()
p.show()
