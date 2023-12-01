import numpy as np
from scipy import constants as cst
import math as m
import matplotlib.pyplot as plt

plt.rcParams.update(plt.rcParamsDefault)
plt.close('all')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": "Computer Modern Roman",
    'font.size': 8
})

colors_table = np.array([['Defect','#53a548'],['Occupancy','#44546A'],['kn','#E6A300'],['kp','#006591'],['en','#FDBD12'],['ep','#6F88FF'],['RSRH','#53a548'],['EF','#C00000'],['EFp','#006591'],['EFn','#E6A300']])
linestyle = [':','-.','-']

                ###########################################
                ### Shockley Read Hall statistics (SRH) ###
                ###########################################

#Flag
savefig_flag = False

#Figure parameters
fontlegend_size = 7
linewidthb5 = 1
linewidthb5_data = 1
markersize_data = 7
labelpadvalue = 3.75
figsize = 5.55226
ticklength = 2.5
alphavalue = 0.05
alphavalueline = 0.25

def kn_func(ET):
    #out = m.exp(1/(ET+0.1))
    out = 1E-20 + m.exp(-1/(0.2*ET))
    return out 

def kp_func(ET):
    #out = m.exp(-1/(ET-1.6))
    out = 1E-20 + m.exp(1/(0.2*(ET-1.5)))
    return out

def n(Vapp,FermiE):
    out = Nc*m.exp(-(GapE-(FermiE+Vapp/2))/kbTeV)
    return out
def p(Vapp,FermiE):
    out = Nv*m.exp(-(FermiE-Vapp/2)/kbTeV)
    return out
def en(kn,ET):
    out = kn*Nc*m.exp((ET-GapE)/(kbTeV))
    return out
def ep(kp,ET):
    out = kp*Nv*m.exp(-ET/(kbTeV))
    return out

def FT(ET,Vapp,kn,kp,FermiE):
    out = (n(Vapp,FermiE)*kn+ep(kp,ET))/(n(Vapp,FermiE)*kn+p(Vapp,FermiE)*kp+en(kn,ET)+ep(kp,ET)) 
    return out

def eta(ET,Vapp,kn,kp,FermiE):
    out = (kn*kp*(n(Vapp,FermiE)*p(Vapp,FermiE)-n(0,FermiE)*p(0,FermiE)))/(n(Vapp,FermiE)*kn+p(Vapp,FermiE)*kp+en(kn,ET)+ep(kp,ET))
    return out 
     
def tau_n(ET,Vapp,kn,kp,FermiE,NT):
    out =  (n(Vapp,FermiE) - n(0,FermiE))/(eta(ET,Vapp,kn,kp,FermiE)*NT)
    return out

def tau_p(ET,Vapp,kn,kp,FermiE,NT):
    out =  (p(Vapp,FermiE) - p(0,FermiE))/(eta(ET,Vapp,kn,kp,FermiE)*NT)
    return out

def compute(kn_init,kp_init,FermiE,Vapp,NT,linestylecompute):
    ET = np.linspace(0,GapE,1000)
    occupancy = np.zeros(len(ET))
    efficiency = np.zeros(len(ET))
    emisssion = np.zeros((2,len(ET)))
    capture = np.zeros((2,len(ET)))
    lifetimes = np.zeros((3,len(ET)))

    for i in range(len(ET)) :
        kp=kp_init#*kp_func(ET[i]) #m3 s-1
        kn=kn_init#*kn_func(ET[i]) #m3 s-1
        occupancy[i] = FT(ET[i],Vapp,kn,kp,FermiE) #/
        emisssion[0][i] = ep(kp,ET[i]) #s-1
        emisssion[1][i] = en(kn,ET[i]) #s-1
        capture[0][i] = n(Vapp,FermiE)*kn #s-1
        capture[1][i] = p(Vapp,FermiE)*kp #s-1
        efficiency[i] = eta(ET[i],Vapp,kn,kp,FermiE) #/
        lifetimes[0][i] = tau_n(ET[i],Vapp,kn,kp,FermiE,NT) #s
        lifetimes[1][i] = tau_p(ET[i],Vapp,kn,kp,FermiE,NT) #s
        lifetimes[2][i] = 1/(1/lifetimes[0][i] + 1/lifetimes[1][i]) #s

    print("n ",n(Vapp,FermiE))
    print("p ",p(Vapp,FermiE))
    print("n0 ",n(0,FermiE))
    print("p0 ",p(0,FermiE))
    ax_SRH.plot(ET,occupancy,linestyle=linestylecompute,color="#11468F",lw=linewidthb5) #label=voltagevalue
    ax_SRH.fill_between(ET,0,occupancy,color="#11468F",alpha = alphavalue)
    ax_SRH.vlines(FermiE,0,1,color=colors_table[np.argwhere(colors_table[:,0]=='EF')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)#label="$E_F$")
    ax_SRH.vlines(FermiE+Vapp/2,0,1,color=colors_table[np.argwhere(colors_table[:,0]=='EFn')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)#label="$E_{F,n}$")
    ax_SRH.vlines(FermiE-Vapp/2,0,1,color=colors_table[np.argwhere(colors_table[:,0]=='EFp')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)#label="$E_{F,p}$")
    #ax_SRH.set_yscale('log')

    ax_SRH_2.plot(ET,emisssion[0],linestyle=linestylecompute,color='#006591',lw=linewidthb5)#,label="$e_p$")
    ax_SRH_2.plot(ET,emisssion[1],linestyle=linestylecompute,color='#ed7d31',lw=linewidthb5)#,label="$e_n$")
    ax_SRH_2.plot(ET,capture[0],linestyle=linestylecompute,color='#C00000',lw=linewidthb5)#,label="$nkn$")
    ax_SRH_2.plot(ET,capture[1],linestyle=linestylecompute,color='#53a548',lw=linewidthb5)#,label="$pkp$")
    ax_SRH_2.vlines(FermiE,np.amin(emisssion),np.amax(emisssion[0]),colors=colors_table[np.argwhere(colors_table[:,0]=='EF')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)
    ax_SRH_2.vlines(FermiE+Vapp/2,np.amin(emisssion),np.amax(emisssion[0]),color=colors_table[np.argwhere(colors_table[:,0]=='EFn')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)
    ax_SRH_2.vlines(FermiE-Vapp/2,np.amin(emisssion),np.amax(emisssion[0]),color=colors_table[np.argwhere(colors_table[:,0]=='EFp')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)
    ax_SRH_2.set_yscale('log')

    NTeff = NT*efficiency/1E6
    ax_SRH_3.plot(ET,NTeff,linestyle=linestylecompute,color='#D82148',lw=linewidthb5)
    ax_SRH_3.vlines(FermiE+Vapp/2,np.amin(lifetimes[1]),np.amax(NTeff),color=colors_table[np.argwhere(colors_table[:,0]=='EFn')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)
    ax_SRH_3.vlines(FermiE-Vapp/2,np.amin(lifetimes[1]),np.amax(NTeff),color=colors_table[np.argwhere(colors_table[:,0]=='EFp')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)
    ax_SRH_3.vlines(FermiE,np.amin(lifetimes[1]),np.amax(NTeff),color=colors_table[np.argwhere(colors_table[:,0]=='EF')][0][0][1],lw=linewidthb5, linestyle = linestylecompute,alpha=alphavalueline)
    #ax_SRH_3.fill_between(ET,np.amin(NTeff),NTeff,color='#D82148',alpha = alphavalue)
    ax_SRH_3.set_yscale('log')

    ax_SRH_4.plot(ET,lifetimes[0],linestyle=linestylecompute,color='#ed7d31',lw=linewidthb5) #electron
    ax_SRH_4.plot(ET,lifetimes[1],linestyle=linestylecompute,color='#006591',lw=linewidthb5) #hole
    #ax_SRH_4.plot(ET,lifetimes[2],linestyle='--',color='white',lw=linewidthb5) #total
    ax_SRH_4.set_yscale('log')

    for a in axes:
        a.tick_params(which='both',length=ticklength,width=linewidthb5,direction='in',pad=5,top=True,right=True)
        a.legend(loc='upper right',frameon=False)
    
    #ax_SRH.set_ylabel("$f_T$",labelpad=4*labelpadvalue)
    ax_SRH.set_ylim(0,1.1)
    ax_SRH.set_xticks([0,0.5,1,1.5])
    #ax_SRH_2.set_ylabel("$nk_n$, $pk_p$, $e_n$, $e_p$ [s$^{-1}$]", labelpad=labelpadvalue)
    #ax_SRH_2.set_ylim(emisssion[1][500],emisssion[1][700])
    #ax_SRH_3.set_ylabel("$R_{RSH}$ [cm$^{-3}$]", labelpad=labelpadvalue)
    #ax_SRH_3.set_ylim(1E5,1E14)
    ax_SRH_3.set_xlabel("$E_T$ [eV]", labelpad=labelpadvalue)
    ax_SRH_3.tick_params(which='both',length=ticklength,width=linewidthb5,direction='in',pad=5,top=True,right=False)
    ax_SRH_3.tick_params()

# Figure generation :
SRHfig, (ax_SRH,ax_SRH_2,ax_SRH_3) = plt.subplots(3,1,figsize=(figsize/3,figsize),sharex = True)
SRHfig.subplots_adjust(left=0.2,bottom=0.075,right=0.8,top=0.995,hspace=0.125)
ax_SRH_4 = ax_SRH_3.twinx()
axes = [ax_SRH,ax_SRH_2,ax_SRH_3,ax_SRH_4]

#Generate graphs
Temperature = 300
GapE = 1.5
kbTeV = cst.Boltzmann*Temperature/cst.elementary_charge
NT = 1E21 #m-3
Nc = 2*m.pow((2*cst.pi*cst.electron_mass*cst.Boltzmann*Temperature)/(cst.Planck*cst.Planck),3/2)
Nv = 2*m.pow((2*cst.pi*cst.electron_mass*cst.Boltzmann*Temperature)/(cst.Planck*cst.Planck),3/2)

Vapvalues = [0.0001,0.5,1]#np.linspace(0.00001,1,3)
print("Voltage ",Vapvalues)
EFvalues = [0.35,0.75,1.05]#np.linspace(0.5,1.2,3)
print("Femi ",EFvalues)
deltakpn = [[2,1,0.5],[0.5,1,2]]
print("Delta ",deltakpn)


kn = 1E-15 #m3 s-1
kp = 1E-15 #m3 s-1
for v,value in enumerate(EFvalues):
    #kp=deltakpn[0][v]*1E-15
    #kn=deltakpn[1][v]*1E-15
    compute(kn,kp,value,0.25,NT,linestyle[v]) 

plt.show()
if savefig_flag:
    SRHfig.savefig('/Users/thomasratz/Desktop/SRH.pdf',dpi=300,facecolor='none',bbox_inches='tight')