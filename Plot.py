import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

data = np.genfromtxt("tmp/CVA.dat",delimiter="\t")
data1 =	np.genfromtxt("tmp/CVAMean.dat",delimiter="\t")
data2 =	np.genfromtxt("tmp/CVAPMean.dat",delimiter="\t")

BF1 =	np.genfromtxt("tmp/BucketFit1.dat",delimiter="\t")
BF2 =	np.genfromtxt("tmp/BucketFit2.dat",delimiter="\t")
BF3 =	np.genfromtxt("tmp/BucketFit3.dat",delimiter="\t")
BF4 =	np.genfromtxt("tmp/BucketFit4.dat",delimiter="\t")
BF5 =	np.genfromtxt("tmp/BucketFit5.dat",delimiter="\t")
BF6 =	np.genfromtxt("tmp/BucketFit6.dat",delimiter="\t")
BF7 =	np.genfromtxt("tmp/BucketFit7.dat",delimiter="\t")
BF8 =	np.genfromtxt("tmp/BucketFit8.dat",delimiter="\t")
BF9 =	np.genfromtxt("tmp/BucketFit9.dat",delimiter="\t")
BF10 =	np.genfromtxt("tmp/BucketFit10.dat",delimiter="\t")
BF11 =	np.genfromtxt("tmp/BucketFit11.dat",delimiter="\t")
BF12 =	np.genfromtxt("tmp/BucketFit12.dat",delimiter="\t")

x = data[:,0]
y = data[:,1]

xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

Mean = np.mean(data1[:,1])

plt.subplot(2,1,1)
plt.scatter(x,y,c=z,s=5, edgecolor='')
# plt.scatter(x,y)
plt.plot(data1[:,0],data1[:,1])
plt.plot(data2[:,0],data2[:,1])
plt.xlabel('Time')
plt.ylabel('Exotic Value ($)')
plt.title('Risk Profile of Exotic Derivative')
plt.annotate(Mean,xy=(1,data[0,1]+0.2*data[0,1]))
plt.grid(True)

plt.subplot(2,1,2)
plt.scatter(BF1[:,1],BF1[:,2])
plt.scatter(BF2[:,1],BF2[:,2])
plt.scatter(BF3[:,1],BF3[:,2])
plt.scatter(BF4[:,1],BF4[:,2])
plt.scatter(BF5[:,1],BF5[:,2])
plt.scatter(BF6[:,1],BF6[:,2])
plt.scatter(BF7[:,1],BF7[:,2])
plt.scatter(BF8[:,1],BF8[:,2])
plt.scatter(BF9[:,1],BF9[:,2])
plt.scatter(BF10[:,1],BF10[:,2])
plt.scatter(BF11[:,1],BF11[:,2])
plt.scatter(BF12[:,1],BF12[:,2])
plt.xlabel('FX Rate')
plt.ylabel('Exotic Value ($)')
plt.title('Value')
plt.grid(True)

plt.show()
