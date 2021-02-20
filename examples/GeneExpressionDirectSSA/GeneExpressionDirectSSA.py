import math, sys, os, time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import bayes_mvs as bayesest

from PyEcoLib.simulator import Simulator


mean_size = 3
doubling_time = 18
tmax = 180
sample_time = 2
div_steps = 10
ncells = 1000

gr = np.log(2)/doubling_time
kd = div_steps * gr/mean_size

sampling_time = sample_time
rprom = 10
pprom = 1000
gammar = 5 * gr
kr = rprom*(gr+gammar)
kp = pprom*gr/rprom

pop = np.zeros([ncells, 6])

indexes = np.int(tmax/sampling_time)
rarray = np.zeros([ncells, indexes])
parray = np.zeros([ncells, indexes])
tarray = np.zeros([indexes])
szarray = np.zeros([ncells, indexes])
cellindex = 0
indexref = 0
start = time.time()
time2 = []
for cell in pop:
    if ncells > 100:
        if cellindex/ncells > indexref:
            print(str(np.int(100*cellindex/ncells))+"%")
            indexref += 0.1

    #Initialize the simulator
    sim = Simulator(ncells=1, gr = gr, sb = mean_size, steps = div_steps)

    #________________________
    #Example of a direct SSA simulation
    cell[0] = mean_size #Initial size
    cell[1] = mean_size*rprom #Initial RNA number
    cell[2] = mean_size*pprom #Initial Protein number
    cell[3] = (1/gr)*np.log(1-(gr/(kr*cell[0]))*np.log(np.random.rand())) #Time to the next RNA creation
    cell[4] = -np.log(np.random.rand())/(gammar*cell[1]) # time to the next RNA degradation
    cell[5] = -np.log(np.random.rand())/(kp*cell[1]) # time to the next protein creation
    t = 0
    reactions=[[0,1,0,0,0,0],[0,-1,0,0,0,0],[0,0,1,0,0,0]] #Reactions (RNA creation, RNA active degradation, Protein creation)
    nextt = 0
    index = 0
    ndiv = 0
    stp = 0
    while t<tmax: #iterating over time
        nr = cell[1]
        nprot = cell[2]
        sz = cell[0]
        tt=sim.get_next_t(0)
        tnextarr = [cell[3],cell[4],cell[5],tt]
        tau = np.min(tnextarr)
        dp = sim.get_dp(0)
        #------------------
        sim.simulate(tmax=tau) #Simulate size dynamics for that given time
        #--------------------
        if np.argmin(tnextarr) != 3:
            cell += reactions[np.argmin(tnextarr)] #if reaction is not a division step

        if sim.get_ndiv(0) > ndiv:
            cell[1] = np.random.binomial(nr,dp) # RNA segregated binomially
            cell[2] = np.random.binomial(nprot,dp) # Protein segregated binomially
            ndiv+=1#stp=0
        cell[0] = sim.get_sz(0)
        nr = cell[1] #Refreshing RNA number
        nprot = cell[2] #Refreshing Protein number
        sz = cell[0] #Refreshing size number
        cell[3] = (1/gr)*np.log(1-(gr/(kr*cell[0]))*np.log(np.random.rand())) #time to thenext rna creation
        cell[4] = -np.log(np.random.rand())/(gammar*cell[1]) #time to the next rna degradation
        cell[5] = -np.log(np.random.rand())/(kp*cell[1]) #time to next protein creation
        t+=tau
        if t > nextt and index<len(tarray): #storing data
            rarray[cellindex,index] = nr/sz # RNA concentration
            parray[cellindex,index] = nprot/sz # Protein concentration
            szarray[cellindex,index] = sz # Cell size
            tarray[index] = t # Time
            index += 1
            nextt += sampling_time
    cellindex += 1
print('It took', np.int(time.time()-start), 'seconds.')


data=pd.DataFrame(np.transpose(np.array(szarray)))
ind=0
newcol=[]
for name in data.columns:
    newcol.append("mom"+str(ind))
    ind+=1
data.columns=newcol
mnszarray=[]
cvszarray=[]
errcv2sz=[]
errmnsz=[]
for m in range(len(data)):
    szs=data.loc[m, :].values.tolist()
    mean_cntr, var_cntr, std_cntr = bayesest(szs,alpha=0.95)
    mnszarray.append(mean_cntr[0])
    errmnsz.append(mean_cntr[1][1]-mean_cntr[0])
    cvszarray.append(var_cntr[0]/mean_cntr[0]**2)
    errv=(var_cntr[1][1]-var_cntr[0])/mean_cntr[0]**2+2*(mean_cntr[1][1]-mean_cntr[0])*var_cntr[0]/mean_cntr[0]**3
    errcv2sz.append(errv)

data['time'] = tarray
data['Mean_sz'] = mnszarray
data['Error_mean'] = errmnsz
data['sz_CV2'] = cvszarray
data['Error_CV2'] = errcv2sz
if not os.path.exists('./data'):
    os.makedirs('./data')
data.to_csv("./data/szsim.csv")

tmax=9*doubling_time
dt=0.0001*doubling_time
lamb=1
a=gr
nsteps=div_steps
k=kd

v0=mean_size
#psz1=[]
ndivs=10
t=0
bigdeltat=0.1
steps=int(np.floor(tmax/dt))
u=np.zeros([ndivs,nsteps])#(DIVS,STEPS)
u[0]=np.zeros(nsteps)
u[0][0]=1#P_00
allmeandivs4=[]#average divisions along the time
allvardiv4=[] # variace of pn along the time
allmeansz4=[]
allvarsz4=[]
time4=[]#time array
yenvol=[]
xenvol=[]
start=0
count=int(np.floor(tmax/(dt*1000)))-1
count2=0
start = time.time()
for l in range(steps):
    utemp=u
    for n in range(len(utemp)):#n=divs,
        for m in range(len(utemp[n])):#m=steps
            if (m==0):#m=steps
                if(n==0):#n=divs
                    dun=-k*v0**lamb*np.exp(lamb*a*t)*(utemp[0][0])
                    u[n][m]+=dun*dt
                else:
                    arg=lamb*(a*t-n*np.log(2))
                    dun=k*v0**lamb*np.exp(arg)*((2**lamb)*utemp[n-1][len(utemp[n])-1]-utemp[n][0])
                    u[n][m]+=dun*dt
            elif(m==len(utemp[n])-1):
                if(n==len(utemp)-1):
                    arg=lamb*(a*t-n*np.log(2))
                    dun=k*v0**lamb*np.exp(arg)*(utemp[n][len(utemp[n])-2])
                    u[n][m]+=dun*dt
                else:
                    arg=lamb*(a*t-n*np.log(2))
                    dun=k*v0**lamb*np.exp(arg)*(utemp[n][m-1]-utemp[n][m])
                    u[n][m]+=dun*dt
            else:
                arg=lamb*(a*t-n*np.log(2))
                dun=k*v0**lamb*np.exp(arg)*(utemp[n][m-1]-utemp[n][m])
                u[n][m]+=dun*dt
    t+=dt
    count=count+1
    if count==int(np.floor(tmax/(dt*1000))):
        time4.append(t/doubling_time)
        mean=0
        for n in range(len(utemp)):
            pnval=np.sum(u[n])
            mean+=n*pnval
        allmeandivs4.append(mean/mean_size)
        var=0
        for n in range(len(utemp)):#divs
            pnval=np.sum(u[n])
            var+=(n-mean)**2*pnval
        allvardiv4.append(np.sqrt(var))
        pn=np.zeros(ndivs)
        sizen=np.zeros(ndivs)
        meansz=0
        for ll in range(len(utemp)):
            pnltemp=np.sum(u[ll])#prob of n divs
            pn[ll]=pnltemp#
            sizen[ll]=np.exp(a*t)/2**ll#
            meansz+=pnltemp*v0*np.exp(a*t)/2**ll
        allmeansz4.append(meansz)
        varsz=0
        for ll in range(len(utemp)):
            pnltemp=np.sum(u[ll])
            varsz+=(v0*np.exp(a*t)/2**ll-meansz)**2*pnltemp
        allvarsz4.append(varsz)
        count=0
        count2+=1
        if(count2==100):
            print(str(int(100*t/tmax))+"%")
            count2=0
print('It took', np.int(time.time()-start), 'seconds.')

fig, ax = plt.subplots(1,2, figsize=(12,4))
#ax[0].plot(tarray,mnszarray)
ax[0].fill_between(np.array(tarray)/doubling_time,np.array(mnszarray)-np.array(errmnsz),np.array(mnszarray)+np.array(errmnsz),
                 alpha=1, edgecolor='#4db8ff', facecolor='#4db8ff',linewidth=0,label='SSA')
#ax[1].plot(tarray,cvszarray)
ax[1].fill_between(np.array(tarray)/doubling_time,np.array(cvszarray)-np.array(errcv2sz),np.array(cvszarray)+np.array(errcv2sz),
                 alpha=1, edgecolor='#4db8ff', facecolor='#4db8ff',linewidth=0)
ax[0].plot(np.array(time4),np.array(allmeansz4),lw=2,c='#006599',label="Numerical")
ax[1].plot(np.array(time4),np.array(allvarsz4)/np.array(allmeansz4)**2,lw=2,c='#006599')

ax[0].set_ylabel("$s$ ($\mu$m)",size=20)
ax[1].set_ylabel("$C_V^2(s)$",size=20)
ax[0].set_xlabel(r"$t/\tau$",size=20)
ax[1].set_xlabel(r"$t/\tau$",size=20)

ax[0].set_ylim([1,1.2*np.max(mnszarray)])
ax[1].set_ylim([0,1.2*np.max(cvszarray)])
for l in [0,1]:
    ax[l].set_xlim([0,tmax/doubling_time])
    taqui=np.arange(0,(tmax+1)/doubling_time,step=1)
    ax[l].set_xticks(np.array(taqui))
    ax[l].grid()
    ax[l].tick_params(axis='x', labelsize=15)
    ax[l].tick_params(axis='y', labelsize=15)
    for axis in ['bottom','left']:
        ax[l].spines[axis].set_linewidth(2)
        ax[l].tick_params(axis='both', width=2,length=6)
    for axis in ['top','right']:
        ax[l].spines[axis].set_linewidth(0)
        ax[l].tick_params(axis='both', width=0,length=6)
plt.subplots_adjust(hspace=0.3,wspace=0.3)
taqui=np.arange(0,0.15,step=0.02)
ax[1].set_yticks(np.array(taqui))
ax[0].legend(fontsize=15)
if not os.path.exists('./figures'):
    os.makedirs('./figures')
plt.savefig('./figures/size_statistics.svg',bbox_inches='tight')
plt.savefig('./figures/size_statistics.png',bbox_inches='tight')
