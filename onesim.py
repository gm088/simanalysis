
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

print("####THIS IS FOR THE ANALYSIS OF ONE distances_tracked.csv FILE and its distances####")
title=input("enter plot title\n")

class interaction:                                       #store information about each tracked interaction
	def __init__(self, line):
	    self.atom = {"index": int(line[25:29]), "type": line[38:41], "resid": int(line[33:38])}       
	    self.ligatom = {"index": int(line[5:9]), "type": line[18:21]}
	    self.ID = line[:2]                               #dunno about this tbh, could be different for same bond in diff replica
	    self.whichtrajs = []
	    self.whichtrajsnot = []

bonds=[]
thing=False
dist_dir="distances"   #directory with distances data

def read_data(distdir):
    CVtrack1=[]
    
    numtrajs = int(os.popen("ls {0}/ | wc -l".format(distdir)).read())-1
    str = (os.popen("wc -l {0}/string_2_distances.dat".format(distdir)).read())
    stride = int(str[0:5])
    
    for i in range(numtrajs):
        traj = np.loadtxt(os.path.join("{0}".format(distdir), "string_{0:d}_distances.dat".format(i+1)), float)
        CVtrack1.append(traj)
    
    CVtrack = np.zeros((numtrajs*traj.shape[0], traj.shape[1]))
    ncols = traj.shape[1]
    
    k=0
    for i in range(numtrajs):
        CVtrack[k:k+CVtrack1[i].shape[0],:] = CVtrack1[i][:,:]
        k += CVtrack1[i].shape[0]
    
    CVtrack1 = None    #free this memory
        
    return CVtrack, stride

dists, strider = read_data(dist_dir)
#the distances we are interested in; find indices from distances_tracked.csv
iwant=[3,7,8,13,20]

################STORING ALL INTERACTION DATA IN A LIST OF CLASS OBJS############

with open("distances_tracked.csv") as f:
    for line in f:
    	if thing:
    		bonds.append(interaction(line))
    	elif line[0]=="I":
    		thing = True
thing = False

##################WHICH PARTs OF THE TRAJ IS THE INTERACTION FROM###############

for i in range(len(bonds)):
    with open("distances_tracked.csv") as f:
    	line0 = f.readline()
    	for j in range(len(line0)):
    		if bonds[i].ID == line0[j-1]+line0[j]:      #big NB - this will have to be changed if bond ID is 3 or more digits
    			index = j                               #initialise index that points to interaction ID
    	for line in f:
    		try:
    			#print(line[index])
    			if line[index]=="X":
    				bonds[i].whichtrajs.append(int(line[5:7]))
    		except IndexError:
    			break

#for i in range(len(bonds)):                                         #trajs in which it is monitored but not biased
#    for j in range(min(bonds[i].whichtrajs),max(bonds[i].whichtrajs)+1,1):
#        if j not in bonds[i].whichtrajs:
#            bonds[i].whichtrajsnot.append(j)

#########NOW WHERE DO WE WRITE TO??############################################

colors = [colors.to_rgba(c)
          for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
colors=colors*5

max_x = max(bonds[0].whichtrajs)
for i in range(len(bonds)):
    if max(bonds[i].whichtrajs) > max_x:
        max_x = max(bonds[i].whichtrajs)
max_x = max_x+2

max_x=24

fig, axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [2.5,1]}, figsize=[10,5])
fig.suptitle(title)

#iwant=range(len(bonds))
kk=0
for i in iwant:
    ind=[]
    for k in range(len(bonds[i].whichtrajs)):
        if bonds[i].whichtrajs[k] == bonds[i].whichtrajs[-1] or bonds[i].whichtrajs[k]+1 != bonds[i].whichtrajs[k+1]:
            if not ind:
                x = np.array(bonds[i].whichtrajs[k])
                x = np.append(x, bonds[i].whichtrajs[k]+0.5)           #I keep forgetting these methods generate a copy of the array
                y = np.full(np.size(x), bonds[i].ID)
                axs[1].plot(x, y, linewidth=1.5, color=colors[kk])
            else:
                ind.append(bonds[i].whichtrajs[k])
                x = np.copy(ind)
                y = np.full(np.size(x), bonds[i].ID)
                axs[1].plot(x, y, linewidth=1.5, color=colors[kk])
                ind=[]
        else:
            ind.append(bonds[i].whichtrajs[k])
    kk+=1

yticks = []
for i in iwant:
    yticks.append("({2:d}){0:s}-{1:s}".format(bonds[i].atom["type"], bonds[i].ligatom["type"], bonds[i].atom["resid"])) 
xticks = []
for i in np.arange(1,max_x):
    xticks.append(range(max_x)[i]*10)   
plt.xlim(0, max_x)
plt.xticks(np.arange(1,max_x), labels=xticks, fontsize=8, rotation=90)
plt.xlabel('time (ns)')
plt.sca(axs[1])
axs[1].set_ylim([-1,len(iwant)])
plt.yticks(range(len(iwant)), labels=yticks, fontsize=5)

for ax in axs:
    ax.label_outer()

#plt.yticks(range(len(bonds)), labels=yticks, fontsize=4.5)
kk=0
dists = dists[0::strider,:]
for i in iwant:
    axs[0].plot(dists[:,i], color=colors[kk])
    kk+=1
axs[0].set_ylabel("CV value")

floorx = np.arange(max_x+1)
floory = np.full(np.size(floorx),3.5)
axs[0].plot(floorx,floory, color='black', linestyle='--')
ceily = np.full(np.size(floorx),7)
axs[0].plot(floorx,ceily, color='black', linestyle='--')
axs[0].set_ylim([1,20])
#axs[1].fill_between(floorx, floory,ceily,alpha=0.2)

plt.savefig("traj_single.png", dpi=250)
