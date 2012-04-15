#
#  this is a simple test bed for playing with FFTs in Python
#
import numpy as np

def filldata(dim1, dim2, data):
   '''fill a data array to use in FFT experiments''' 
   import math
   dc = 0.75
   f1 = 2.0
   f2 = 7.0
   for x in range(dim1):
       for y in range(dim2):
          data[x,y] = dc + math.cos(x)*math.sin(y)

def compare(dim1, dim2, data1, data2):
   '''compute sum of differences squared between two arrays'''
   errsq = 0.0
   for x in range(dim1):
       for y in range(dim2):
          errsq = errsq + (data1[x,y]-data2[x,y])**2
   return errsq

d1 = 4
d2 = 4
d = np.zeros([d1,d2], dtype=np.float)
#f = np.zeros([d1,d2], dtype=np.float)
g = np.zeros([d1,d2], dtype=np.float)

filldata(d1, d2, d)

f = np.fft.fftn(d)
g = np.real_if_close(np.fft.ifftn(f), tol = 100000)

print d

print g



errsq = compare(d1, d1, d, g)

print("error squared from FFTs is ",errsq)

