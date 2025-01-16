import numpy as np
import sys
import scipy.constants
import matplotlib.pyplot as plt

Na = scipy.constants.Avogadro # Avogadro's number
kb = scipy.constants.k # Boltzmann Constant k
data1 = sys.argv[1]
data2 = sys.argv[2]
run = sys.argv[3]
temp = sys.argv[4]
#data3 = sys.argv[3]
#data4 = sys.argv[4]
#data5 = sys.argv[5]
#data6 = sys.argv[6]
#run = sys.argv[7]
#temp = sys.argv[8]
dataa = []
data1a = np.loadtxt(fname = data1,skiprows = 17,usecols = 1)
data2a = np.loadtxt(fname = data2,skiprows = 17,usecols = 1)
#data3a = np.loadtxt(fname = data3,skiprows = 17,usecols = 1)
#data4a = np.loadtxt(fname = data4,skiprows = 17,usecols = 1)
#data5a = np.loadtxt(fname = data5,skiprows = 17,usecols = 1)
#data6a = np.loadtxt(fname = data6,skiprows = 17,usecols = 1)
#dataa = (data1a + data2a + data3a + data4a + data5a + data6a)/6
dataa = (data1a + data2a)/2
out = open("vac-CH4-"+temp+"K-"+run+"-average-prod1.txt","w")
#out = open("vac-Ar-"+temp+"K-"+run+"-average-prod1.txt","w")
out.truncate()
for i in range(len(dataa)):
	print(float(dataa[i]), file = out)
out.close()
print("done!")
