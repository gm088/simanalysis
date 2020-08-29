from pyemma.thermo import tram
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

numsims = int(os.popen('ls data/ | wc -l').read())

data = []

for i in range(numsims):
	traj = np.loadtxt(os.path.join("data", "data_s100_string_{0:d}.dat".format(i+1)), dtype=np.double)
    #stride if u want
	traj = traj[0:10000:10,:]
	data.append(traj)

potcentres = np.loadtxt(os.path.join("..","optstring","newconstr_3.dat"), dtype=np.double)
k_val = 20
datlength = traj.shape[0]
#datlength = 10000
numsims = len(data)

q = np.zeros((numsims,datlength), dtype=np.double)

for i in range(numsims):
    for j in range(datlength):
        q[i,j] += sum(data[i][j,1:8])
        
bias = []

for i in range(numsims):
    abias = []
    for j in range(datlength):
        ubias = []
        for k in range(numsims):
            UB = sum(0.5*k_val*((data[i][j,:] - potcentres[k,:])**2))
            ubias.append(UB)
        abias.append(ubias)
    bias.append(abias)
    
for i in range(200):
    bias[i] = np.asarray(bias[i], dtype=np.double)


v_max = np.amax(q)
v_min = np.amin(q)
ttrajs = []
dtrajs = []

edges = np.linspace(v_min, v_max, 200)
for i in range(numsims):
    inds = np.digitize(q[i,:], edges)
    dtrajs.append(inds)
for i in range(numsims):
    ttrajs.append(np.full((datlength),i))


tram_obj = tram(ttrajs, dtrajs, bias, 1)

tram_obj.save('tramstr3')

plt.plot(tram_obj.f)
plt.savefig("fren.png", dpi=199)
