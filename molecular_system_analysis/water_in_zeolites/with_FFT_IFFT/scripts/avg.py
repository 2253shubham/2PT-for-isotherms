import numpy as np
import sys

data1 = sys.argv[1]
data2 = sys.argv[2]
data4 = sys.argv[3]
data3 = []
data1a = np.loadtxt(fname = data1)
data1a = data1a.reshape(len(data1a),-1)
data2a = np.loadtxt(fname = data2)
data2a = data1a.reshape(len(data2a),-1)
data3 = (data1a+data2a)/2
out = open(data4,"w")
out.truncate()
for i in range(len(data3)):
    print(*data3[i], file = out)
out.close()


