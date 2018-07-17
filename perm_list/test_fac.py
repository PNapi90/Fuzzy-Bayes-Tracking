from numpy import *

def fac(a):
    retval = 1
    for i in range(1,a+1):
        retval *= i
    return retval




a = 0
for i in range(1,11):
    a += fac(i)*i

print(a*4/(1024*1024))