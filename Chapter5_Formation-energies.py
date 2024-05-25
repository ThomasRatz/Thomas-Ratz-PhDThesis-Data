                                ###################################################
                                #                                                 #
                                #       Formation energy of defects VS E_F        #   
                                #                                                 #
                                ###################################################

# Librairies - Python version 
import string
import numpy as np
import matplotlib.pyplot as plt # Enable the plot function
from matplotlib.widgets import Slider
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import math as m
from scipy import constants as cst
#from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
import seaborn as sns
font = {'family': 'serif', 'size': 18}
fontlabel = {'family': 'serif', 'size': 16}
fontlabelbar = {'family': 'serif', 'size': 12}
axcolor = 'lightgrey'
plt.rcParams['axes.linewidth'] = 1.5

#Parameters
Set_Slider = True                       #Slider for chemical potential values in formation energy plot
draweachdefect = False          
CZGS_calc = True                       #Compute CZGS point defects
savefigs = True                        #Save figures in defined "pathsave"
pathsave = '/Users/thomasratz/Desktop/'
oneelectrontransitiononly = False       #Draw only transition difference equals to one electron charge
yenergy = 6                             #Formation energy range for \DeltaH_F [eV]
nline = 1000                            #Points per line in \DeltaH_F
Madelung = 0.304                        #Madelung constant for correction term

#CZTS --> 0.46696696696696693/((nbkb*cst.Boltzmann*300)/cst.elementary_charge) EF equilibrium
#CZGS --> 0.39444444444444443/((nbkb*cst.Boltzmann*300)/cst.elementary_charge)
nbkb = 0.46696696696696693/((cst.Boltzmann*300)/cst.elementary_charge)
nbkb_CZGS = 0.39444444444444443/((cst.Boltzmann*300)/cst.elementary_charge)

#point E in phase diagram Fig.5.1 for CZTS
mu_cu = -0.550#walsh[0]#adfm[2][0]#-0.55
mu_zn = -1.558#walsh[1]#adfm[2][1]#-1.56
mu_sn = -0.561#walsh[2]#adfm[2][2]#-0.56
chemicalpotentialpath_CZTS = '/Users/thomasratz/Desktop/path_CZTS.txt'  #Chemical potential path for Fig.5.3

#point E in phase diagram Fig.5.1 for CZGS
CZGS_mu_cu = -0.550 #point E
CZGS_mu_zn = -1.506
CZGS_mu_ge = -0.675
chemicalpotentialpath_CZGS = '/Users/thomasratz/Desktop/path_CZGS.txt' #Chemical potential path for Fig.5.3

########################################################################################################################
#                                                                                                                      #
#                                           Formation energy plot for CZTS                                             #
#                                                                                                                      #
########################################################################################################################
#Calculation data
Pure_energies_CZTS = [["Cu",-14.58743314,-3.65226],["Zn",-2.53611365,-1.2574],["Sn",-36.27147174,-4.53397],["S",-167.63710169,-5.23831]]
Formation_CZTS = 8*(-308.96612677/64) - 2 * Pure_energies_CZTS[0][2] - Pure_energies_CZTS[1][2] - Pure_energies_CZTS[2][2] - 4*Pure_energies_CZTS[3][2]

#[Energy, VBM, Eg]
#CZTS_host_8 = [,,1.316] 
#Fermi in DOSCAR CZTS supercell = 4.12943529 eV and band gap 1.179442
#Fermi in CZTS 8 atoms : 3.863666 and band gap 1.316
CZTS_host = [-308.96612677,4.12943529,1.316] 

#[0-"Name", 1-charge state, 2-n_Cu, 3-n_Zn, 4-n_X, 5-n_S, 6-total energy, 7-Formation energy, 8-Potential correction compare to host supercell]
CZTS_Defect = [[["$V_{Cu}$",0,-1,0,0,0,-303.757204199,[0,0],0.08393984],["$V_{Cu}^{-1}$",-1,-1,0,0,0,-300.062440479,[0,0],0.16932281]],
        [["$V_{Sn}$",0,0,0,-1,0,-298.669605883,[0,0],0.23182653],["$V_{Sn}^{-1}$",-1,0,0,-1,0,-294.887822033,[0,0],0.33732134],["$V_{Sn}^{-2}$",-2,0,0,-1,0,-291.099004956,[0,0],0.44629551],["$V_{Sn}^{-3}$",-3,0,0,-1,0,-287.512231044,[0,0],0.55045844],["$V_{Sn}^{-4}$",-4,0,0,-1,0,-283.933338404,[0,0],0.65582942]],
        [["$V_{Zn}$",0,0,-1,0,0,-304.026156815,[0,0],0.05121121],["$V_{Zn}^{-1}$",-1,0,-1,0,0,-300.250930989,[0,0],0.14345909],["$V_{Zn}^{-2}$",-2,0,-1,0,0,-296.589002191,[0,0],0.23541598]],
        [["$V_S$",0,0,0,0,-1,-301.713872099,[0,0],0.21667656],["$V_S^{+1}$",+1,0,0,0,-1,-305.727824725,[0,0],0.08321605]],
        [["$Cu_{Zn}$",0,+1,-1,0,0,-309.716311444,[0,0],-0.04339986],["$Cu_{Zn}^{-1}$",-1,+1,-1,0,0,-305.632559386,[0,0],0.03499529]],
        [["$Zn_{Cu}$",0,-1,+1,0,0,-306.264246427,[0,0],0.04574733],["$Zn_{Cu}^{+1}$",+1,-1,+1,0,0,-312.011885874,[0,0],-0.00469694]], #312.011885874
        [["$Cu_{Sn}$",0,+1,0,-1,0,-305.183720226,[0,0],0.09792049],["$Cu_{Sn}^{-1}$",-1,+1,0,-1,0,-301.155552755,[0,0],0.18634445],["$Cu_{Sn}^{-2}$",-2,+1,0,-1,0,-297.231915182,[0,0],0.27995114],["$Cu_{Sn}^{-3}$",-3,+1,0,-1,0,-293.207109971,[0,0],0.38878789]],
        [["$Sn_{Cu}$",0,-1,0,+1,0,-307.090186587,[0,0],-0.02829221],["$Sn_{Cu}^{+1}$",+1,-1,0,+1,0,-312.919053527,[0,0],-0.08220638],["$Sn_{Cu}^{+2}$",+2,-1,0,+1,0,-317.601352675,[0,0],-0.09518029],["$Sn_{Cu}^{+3}$",+3,-1,0,+1,0,-322.955792179,[0,0],0.08952915]],
        [["$Zn_{Sn}$",0,0,+1,-1,0,-304.573017074,[0,0],0.15147367],["$Zn_{Sn}^{-1}$",-1,0,+1,-1,0,-300.862483469,[0,0],0.23167036],["$Zn_{Sn}^{-2}$",-2,0,+1,-1,0,-297.157372061,[0,0],0.31188925]],
        [["$Sn_{Zn}$",0,0,-1,+1,0,-309.766179707,[0,0],-0.07180158],["$Sn_{Zn}^{+1}$",1,0,-1,+1,0,-314.700235124,[0,0],-0.12648315],["$Sn_{Zn}^{+2}$",2,0,-1,+1,0,-320.210855316,[0,0],-0.02182494]],
        [["$Cu_{I}$",0,+1,0,0,0,-310.747582263,[0,0],-0.0361141],["$I_{Cu}^{+1}$",+1,+1,0,0,0,-316.146064111,[0,0],-0.00135404]],
        [["$Sn_{I}$",0,0,0,+1,0,-309.439347709,[0,0],-0.05682777],["$I_{Sn}^{+1}$",+1,0,0,+1,0,-314.067843316,[0,0],-0.16879601],["$I_{Sn}^{+2}$",+2,0,0,+1,0,-319.623455527,[0,0],-0.21260955],["$I_{Sn}^{+3}$",+3,0,0,+1,0,-324.270049687,[0,0],-0.26286416],["$I_{Sn}^{+4}$",+4,0,0,+1,0,-329.457289221,[0,0],-2.31646575]],
        [["$Zn_{I}$",0,0,+1,0,0,-308.702426606,[0,0],0.06426412],["$I_{Zn}^{+1}$",+1,0,+1,0,0,-313.077640108,[0,0],-0.01399986],["$I_{Zn}^{+2}$",+2,0,+1,0,0,-317.769505936,[0,0],-2.0719705]],
        [["$S_I$",0,0,0,0,+1,-312.12476178,[0,0],1.13137349],["$I_S^{-1}$",-1,0,0,0,+1,0,[0,0],-0.74563827]]]

#[0-"Name", 1-charge state, 2-n_Cu, 3-n_Zn, 4-n_X, 5-n_S, 6-total energy, 7-Formation energy, 8-Potential correction compare to host supercell]
CZTS_doping = [[["$Ge_{Cu}$",0,-1,0,0,0,-307.734830483,[0,0],0.03547223],["$Ge_{Cu}^{+1}$",+1,-1,0,0,0,-313.527302697,[0,0],-0.01297549],["$Ge_{Cu}^{+2}$",+2,-1,0,0,0,-318.048896958,[0,0],0.04968672],["$Ge_{Cu}^{+1}$",+3,-1,0,0,0,-323.425357744,[0,0],0.25928952]],
               [["$Ge_{Zn}$",0,0,-1,0,0,-310.47592843,[0,0],-0.01075182],["$Ge_{Zn}^{+1}$",+1,0,-1,0,0,-315.274162501,[0,0],0.04632261],["$Ge_{Zn}^{+2}$",+2,0,-1,0,0,-320.889903527,[0,0],0.10314709]],
               [["$Ge_{Sn}$",0,0,0,-1,0,-309.72837172,[0,0],0.15588854]],
               #[["$Ge_{S}$",0,0,0,0,-1,-305.97335478,[0,0],0.01844333]],
               [["$Ge_{I}$",0,0,0,0,0,-310.837914103,[0,0],0.03368119],["$Ge_{I}^{+1}$",1,0,0,0,0,-315.81683,[0,0],-0.02387102],["$Ge_{I}^{+2}$",2,0,0,0,0,-321.501716898,[0,0],-0.0632061],["$Ge_{I}^{+3}$",3,0,0,0,0,-325.904743882,[0,0],-0.14488805],["$Ge_{I}^{+4}$",4,0,0,0,0,-330.372890156,[0,0],-0.22849945]]]

#C_sh, Epsilon, Alpha_M, L
Madelung_CZTS = [-0.342,6.77*cst.epsilon_0,Madelung,10.94234*1E-10] #1.63664

fig_CZTS_all, (ax_CZTS_V, ax_CZTS_S, ax_CZTS_I) = plt.subplots(1,3,figsize=(12,7))
fig_CZTS_all.subplots_adjust(bottom=0.3, wspace = 0.2,right=0.98,left=0.1,top=0.98)
plot_color = ["#00708F","#E6A100","#E61F00","#4D009B","#B98000","#00A053","#00BEB6","#005228","#754D00","#720008","#FF4365","#2B50AA","#FF9FE5","#73A6AD","#FA7921","#CB48B7","#016FB9","#7776BC","#40476D"]
subplotarry = [ax_CZTS_V,ax_CZTS_S,ax_CZTS_I]

#-----------------------------------------------------------------------------------------------------------------------------#
#                                               Formation Energy Calculations                                                 #
#-----------------------------------------------------------------------------------------------------------------------------#

def Ge_doping_potential(mu_Cu,mu_Zn,mu_Sn,mu_S): #calculation of the Ge chemical potential for doping inside Sn-kesterite

        Secondaries_Ge = [["Cu2GeS3",-123.48186136,-123.48186136/24],["Cu3Ge",-32.92180255,-32.92180255/8],["GeS",-44.13581527,-44.13581527/8],
                          ["GeS2",-272.26921339,-272.26921339/48],["SnGeS3",-0.110674289220E+03,-0.110674289220E+03/20]]
        Pure_energies_CZGS = [["Cu",-14.58743314,-3.65226],["Zn",-2.53611365,-1.2574],["Ge",-36.27147174,-5.35416],["S",-167.63710169,-5.23831],["Sn",-36.27147174,-4.53397]]
        Secondaries_formation_Ge= [["Cu2GeS3",6*Secondaries_Ge[0][2] - 2*Pure_energies_CZGS[0][2] -Pure_energies_CZGS[2][2] - 3*Pure_energies_CZGS[3][2] - 0.27],
                                   ["Cu3Ge",4*Secondaries_Ge[1][2] - 3*Pure_energies_CZGS[0][2] - Pure_energies_CZGS[2][2] - 0.27],
                                   ["GeS",2*Secondaries_Ge[2][2] - Pure_energies_CZGS[2][2] - Pure_energies_CZGS[3][2] - 0.27],
                                   ["GeS2",3*Secondaries_Ge[3][2] - Pure_energies_CZGS[2][2] - 2*Pure_energies_CZGS[3][2] - 0.27],
                                   ["SnGeS3",5*Secondaries_Ge[4][2] - Pure_energies_CZGS[4][2] - Pure_energies_CZGS[2][2] - 3*Pure_energies_CZGS[3][2] - 0.27 ]]
        mu_Ge_rich = -0.01
        flag = True
        while flag:
                if(2*mu_Cu + mu_Ge_rich + 3*mu_S <= Secondaries_formation_Ge[0][1]
                  and 3*mu_Cu + mu_Ge_rich <= Secondaries_formation_Ge[1][1]
                  and mu_Ge_rich + mu_S <= Secondaries_formation_Ge[2][1]
                  and mu_Ge_rich + 2*mu_S <= Secondaries_formation_Ge[3][1]
                  and mu_Sn + mu_Ge_rich + 3* mu_S <= Secondaries_formation_Ge[4][1]):
                        flag = False
                mu_Ge_rich -= 0.01
        return mu_Ge_rich

def CZTS_Defect_calc_all(mu_cu,mu_zn,mu_sn,Fermi,x,type,charge,Madelung): # Function to compute the defect formation energy at a given energy "x" (Fermi position)
        CZTS_chemical = [mu_cu,mu_zn,mu_sn,0] #mu_Cu mu_Zn mu_Sn
        CZTS_chemical[3] = (Formation_CZTS - 2*CZTS_chemical[0] - CZTS_chemical[1] - CZTS_chemical[2])/4 #mu_S
        out = 0
        #Hf    #E(i,h)                       #Ehost         #charge                       #VBM         # Madelung electrostatic interaction                                                                                                         #Potentiel alignement   
        out = CZTS_Defect[type][charge][6] - CZTS_host[0] + CZTS_Defect[type][charge][1]*(Fermi + x) + ((0.65)*(m.pow(CZTS_Defect[type][charge][1]*cst.elementary_charge,2)*Madelung)/(2*Madelung_CZTS[1]*Madelung_CZTS[3]))/cst.elementary_charge - CZTS_Defect[type][charge][1]*CZTS_Defect[type][charge][8]
        for j in range(4):
                #Hf          #ni                             #mu_i             #E_i
                out = out - CZTS_Defect[type][charge][2+j]*(CZTS_chemical[j] + Pure_energies_CZTS[j][2])
        return out

def CZTS_Defect_doping(mu_cu,mu_zn,mu_sn,Fermi,x,type,charge,Madelung): # Function to compute the Ge doping defect formation energy at a given energy "x" (Fermi position)
        CZTS_chemical = [mu_cu,mu_zn,mu_sn,0] #mu_Cu mu_Zn mu_Sn
        CZTS_chemical[3] = (Formation_CZTS - 2*CZTS_chemical[0] - CZTS_chemical[1] - CZTS_chemical[2])/4 #mu_S
        Ge_doping_chemical_pot = 0
        Ge_doping_chemical_pot = Ge_doping_potential(CZTS_chemical[0],CZTS_chemical[1],CZTS_chemical[2],CZTS_chemical[3])
        out = 0
        Pure_energies_Ge = -5.35416
        #Hf    #E(i,h)                       #Ehost         #Ge dopant                                    #charge                       #VBM         # Madelung electrostatic interaction                                                                                                         #Potentiel alignement   
        out = CZTS_doping[type][charge][6] - CZTS_host[0] - (Ge_doping_chemical_pot + Pure_energies_Ge) + CZTS_doping[type][charge][1]*(Fermi + x) + ((0.65)*(m.pow(CZTS_doping[type][charge][1]*cst.elementary_charge,2)*Madelung)/(2*Madelung_CZTS[1]*Madelung_CZTS[3]))/cst.elementary_charge - CZTS_doping[type][charge][1]*CZTS_doping[type][charge][8]
        for j in range(4):
                #Hf          #ni                             #mu_i             #E_i
                out = out - CZTS_doping[type][charge][2+j]*(CZTS_chemical[j] + Pure_energies_CZTS[j][2])
        return out

#-----------------------------------------------------------------------------------------------------------------------------#
#                                              Defect Plot Function                                                           #
#-----------------------------------------------------------------------------------------------------------------------------#

def CZTS_Defect_plot(mu_cu,mu_zn,mu_sn,VBM,Madelung): # Plot function for defect formation energy VS Fermi energy
        init = -1 # EF init value (within valence band)
        end = CZTS_host[2]+1  # EF end value (within conduction band)
        x = np.linspace(init,end,nline) #1000 pts from init to end for each line
        findEfMin = np.zeros((len(CZTS_Defect),len(x)))

        for defect,defectvalue in enumerate(CZTS_Defect): # Defect type for intrinsic doping
                # Compute formation energy for one defect type, in each charge state and output the minimal formation energy for each fermi energy value, then plot it
                line_defect_out = np.zeros(len(x))
                line_defect = np.zeros(len(x))
                line_defect_temp = np.zeros(len(x))

                # Formation energy at a given energy E
                for charge, chargevalue in enumerate(CZTS_Defect[defect]): # Defect charge state
                        for a in range(len(x)): # Loop over the energies
                                line_defect_temp[a] = CZTS_Defect_calc_all(mu_cu,mu_zn,mu_sn,VBM,x[a],defect,charge,Madelung) # Calculate the formation energy
                        line_defect = np.vstack((line_defect,line_defect_temp)) # Put 2 lines arrays in a same arrays as [[...,...],[...,...]]
                line_defect = np.delete(line_defect, 0, 0) # Remove first line which is all 0 elements 
                # Find the minimum of all lines (i.e. minimum of formation energy for each charge state of one defect) and also gives the energy at which the deviation occurs
                line_defect_out = np.min(line_defect, axis=0)
                findEfMin[defect] = line_defect_out

                if(defect <= 3): #<=3     
                        ax_CZTS_V.plot(x,line_defect_out,lw = 1, label=CZTS_Defect[defect][0][0], color=plot_color[defect])
                        if draweachdefect == True:
                                for charge in range(len(CZTS_Defect[defect])): # Defect charge state    
                                        ax_CZTS_V.plot(x,line_defect[charge],lw = 0.5, color=plot_color[defect],linestyle='dashed')
                elif(defect > 3 and defect <= 9):    
                        ax_CZTS_S.plot(x,line_defect_out,lw = 1, label=CZTS_Defect[defect][0][0], color=plot_color[defect])
                        if draweachdefect == True:
                                for charge in range(len(CZTS_Defect[defect])): # Defect charge state    
                                        ax_CZTS_S.plot(x,line_defect[charge],lw = 0.5, color=plot_color[defect],linestyle='dashed')
                elif(defect > 9):
                        ax_CZTS_I.plot(x,line_defect_out,lw = 1, label=CZTS_Defect[defect][0][0], color=plot_color[defect])
                        if draweachdefect == True:
                                for charge in range(len(CZTS_Defect[defect])): # Defect charge state    
                                        ax_CZTS_I.plot(x,line_defect[charge],lw = 0.5, color=plot_color[defect],linestyle='dashed')
        
        findEfMin2 = np.zeros(len(x))
        for i,ivalue in enumerate(x):
                for j in range(len(CZTS_Defect)):
                        findEfMin2[i] += findEfMin[j][i]
        print("Defect plot")
        for defect,defectvalue in enumerate(CZTS_doping): # Defect type For Ge extrinsic doping
                # Compute formation energy for one defect type, in each charge state and output the minimal formation energy for each fermi energy value, then plot it
                line_defect_out = np.zeros(len(x))
                line_defect = np.zeros(len(x))
                line_defect_temp = np.zeros(len(x))
                # Formation energy at a given energy E
                for charge,chargevalue in enumerate(CZTS_doping[defect]): # Defect charge state
                        for a in range(len(x)): # Loop over the energies
                                line_defect_temp[a] = CZTS_Defect_doping(mu_cu,mu_zn,mu_sn,VBM,x[a],defect,charge,Madelung) # Calculate the formation energy
                        line_defect = np.vstack((line_defect,line_defect_temp)) # Put 2 lines arrays in a same arrays as [[...,...],[...,...]]
                line_defect = np.delete(line_defect, 0, 0) # Remove first line which is all 0 elements 
                # Find the minimum of all lines (i.e. minimum of formation energy for each charge state of one defect) and also gives the energy at which the deviation occurs
                line_defect_out = np.min(line_defect, axis=0)

                if(defect < 3):    
                        ax_CZTS_S.plot(x,line_defect_out,lw = 1, label=CZTS_doping[defect][0][0], color=plot_color[14+defect])
                        if draweachdefect == True:
                                for charge in range(len(CZTS_doping[defect])): # Defect charge state    
                                        ax_CZTS_S.plot(x,line_defect[charge],lw = 0.5, color=plot_color[14+defect],linestyle='dashed')
                else :
                        ax_CZTS_I.plot(x,line_defect_out,lw = 1, label=CZTS_doping[defect][0][0], color=plot_color[14+defect])
                        if draweachdefect == True:
                                for charge in range(len(CZTS_doping[defect])): # Defect charge state    
                                        ax_CZTS_I.plot(x,line_defect[charge],lw = 0.5, color=plot_color[14+defect],linestyle='dashed')
                              
        for plotdefect in subplotarry:
            plotdefect.vlines(x=0, ymin=-10, ymax=10, color="navy", lw=0.75)                                            # label="$\epsilon_{VBM}$"
            plotdefect.vlines(x=CZTS_host[2], ymin=-10, ymax=10, color="darkorange", lw=0.75)
            plotdefect.vlines(x=(nbkb*cst.Boltzmann*300)/cst.elementary_charge, ymin=-10, ymax=10, color="maroon",lw=1.5, linestyle = '--')   # label="$E_{G}$"
            plotdefect.hlines(y=0, xmin=-10, xmax=10, color="black", lw=0.25,linestyle='dashed')                        # 0 formation energy line
            plotdefect.fill_between(np.arange(-10,0.01,0.01),-10,10,color="navy", alpha=0.1)                            # Conduction band
            plotdefect.fill_between(np.arange(CZTS_host[2],CZTS_host[2]+5,0.01),-10,10,color="darkorange", alpha=0.1)   # Valence band

        ax_CZTS_V.set_ylabel("$\Delta H_F$ [eV]",fontdict = fontlabel , labelpad=5)
        for plotdefect in subplotarry:
            plotdefect.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)
            plotdefect.legend(loc='upper right',prop={'size': 14})
            plotdefect.set_xlabel("$E_F$ [eV]",fontdict = fontlabel , labelpad=5)
            plotdefect.set_xlim(-0.5,CZTS_host[2]+0.5)
            plotdefect.set_ylim(0,yenergy)
            plotdefect.set_yticks(np.arange(-0.5,yenergy+0.5,0.5))
        ax_CZTS_S.legend(loc='upper right',prop={'size': 14},ncol=2,handleheight=2, labelspacing=0.01)
CZTS_Defect_plot(mu_cu,mu_zn,mu_sn,CZTS_host[1],Madelung)    

#-----------------------------------------------------------------------------------------------------------------------------#
#                                                 Ionization  Levels                                                          #
#-----------------------------------------------------------------------------------------------------------------------------#
#The ionization levels is defined as the the Fermi energy for which a defect formation energy of a defect "a->CZTS_Defect[i]"
#in a charge state "q->CZTS_Defect[i][j]" is equal to the defect formation energy of the same defect "a->CZTS_Defect[i]" in 
#another charge state q' -> use intersection function 
defect_number = len(CZTS_Defect)+len(CZTS_doping)
energy_broadening = 0 

fig_CZTS_ionization = plt.figure(figsize=(10,5))
fig_CZTS_ionization.subplots_adjust(left=0.121,bottom=0.138,right=1,top=0.97)
ax_CZTS_ionization = plt.subplot()

CZTS_Defect_label = ["" for x in range(len(CZTS_Defect)+len(CZTS_doping))] # Set defect name
for i in range(len(CZTS_Defect)):
        CZTS_Defect_label[i] = CZTS_Defect[i][0][0]
for i in range(len(CZTS_doping)):
        CZTS_Defect_label[14+i] = CZTS_doping[i][0][0]

def line_intersection(line1,line2): #Find intersection point between two different charge states of a defect (between two lines)
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       return CZTS_host[2]+0.1,0 #Check this !

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

#Function to look between same defects the various charge states transition energies -> use function lineintersection
#between charge 1 and 2,3,4,5,...N
#between charge 2 and 3,4,5 N
ionisationenergies = []
def CZTS_transition(mu_cu,mu_zn,mu_sn,VBM,Madelung):
        color_ionisation = sns.color_palette("magma", as_cmap=True)
        CZTS_chemical = [mu_cu,mu_zn,mu_sn,0] #mu_Cu mu_Zn mu_Sn mu_S
        CZTS_chemical[3] = (Formation_CZTS - 2*CZTS_chemical[0] - CZTS_chemical[1] - CZTS_chemical[2])/4 #mu_S
        for type,typevalue in enumerate(CZTS_Defect): # Defect type
                for charge in range(len(CZTS_Defect[type])): # Defect charge state
                        CZTS_Defect[type][charge][7][0] = CZTS_Defect_calc_all(mu_cu,mu_zn,mu_sn,VBM,0,type,charge,Madelung)
                        CZTS_Defect[type][charge][7][1] = CZTS_Defect_calc_all(mu_cu,mu_zn,mu_sn,VBM,CZTS_host[2],type,charge,Madelung)
        for type,typevalue in enumerate(CZTS_Defect):
                ionisationenergies.append(CZTS_Defect[type][0][0])
                #Search for intersection between defects formation energies:
                for i,ivalue in enumerate(CZTS_Defect[type]):# Defects types
                        for j,jvalue in enumerate(CZTS_Defect[type]): # Defects charges
                                if oneelectrontransitiononly == True:
                                        if i != j and i != i+j and abs(i-j) == 1: #last argument is to take into account only 1e transition
                                                ionisationenergies.append([i,i+j])
                                                #look for interesction between two consecutive charged defects
                                                lines_intersect = [[[-energy_broadening,CZTS_Defect[type][i][7][0]],[CZTS_host[2]+energy_broadening,CZTS_Defect[type][i][7][1]]],
                                                                [[-energy_broadening,CZTS_Defect[type][j][7][0]],[CZTS_host[2]+energy_broadening,CZTS_Defect[type][j][7][1]]]]
                                                #look for intersection and write in table ["type",[i,i+j+1],x_intersection,y_intersection]
                                                ionisationenergies.append(line_intersection(lines_intersect[0],lines_intersect[1]))
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],type-0.41,type+0.41,color='black',lw=3.5)
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],type-0.4,type+0.4,color=color_ionisation(line_intersection(lines_intersect[0],lines_intersect[1])[1]/4),lw=3)
                                                if type <= 3: #<=3
                                                        ax_CZTS_V.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'o', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                elif type > 3 and type <= 9 :
                                                        ax_CZTS_S.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'s', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                elif type > 9 :
                                                        ax_CZTS_I.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'D', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                else :
                                        if i != j and i != i+j :
                                                ionisationenergies.append([i,i+j])
                                                #look for interesction between two consecutive charged defects
                                                lines_intersect = [[[-energy_broadening,CZTS_Defect[type][i][7][0]],[CZTS_host[2]+energy_broadening,CZTS_Defect[type][i][7][1]]],
                                                                [[-energy_broadening,CZTS_Defect[type][j][7][0]],[CZTS_host[2]+energy_broadening,CZTS_Defect[type][j][7][1]]]]
                                                #look for intersection and write in table ["type",[i,i+j+1],x_intersection,y_intersection]
                                                ionisationenergies.append(line_intersection(lines_intersect[0],lines_intersect[1]))
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],type-0.41,type+0.41,color='black',lw=3.5)
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],type-0.4,type+0.4,color=color_ionisation(line_intersection(lines_intersect[0],lines_intersect[1])[1]/4),lw=3)
                                                if type <= 3:#<=3
                                                        ax_CZTS_V.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'o', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                elif type > 3 and type <= 9 :
                                                        ax_CZTS_S.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'s', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                elif type > 9 :
                                                        ax_CZTS_I.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'D', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
        
        print("ionisation plot")
        for type,typevalue in enumerate(CZTS_doping): # Defect type
                for charge in range(len(CZTS_doping[type])): # Defect charge state
                        CZTS_doping[type][charge][7][0] = CZTS_Defect_doping(mu_cu,mu_zn,mu_sn,VBM,0,type,charge,Madelung)
                        CZTS_doping[type][charge][7][1] = CZTS_Defect_doping(mu_cu,mu_zn,mu_sn,VBM,CZTS_host[2],type,charge,Madelung)
        for type,typevalue in enumerate(CZTS_doping):
                ionisationenergies.append(CZTS_doping[type][0][0])
                #Search for intersection between defects formation energies:
                for i,ivalue in enumerate(CZTS_doping[type]):# Defects types
                        for j,jvalue in enumerate(CZTS_doping[type]): # Defects charges
                                if oneelectrontransitiononly == True:
                                        if i != j and i != i+j and abs(i-j) == 1: #last argument is to take into account only 1e transition
                                                ionisationenergies.append([i,i+j])
                                                #look for interesction between two consecutive charged defects
                                                lines_intersect = [[[-energy_broadening,CZTS_doping[type][i][7][0]],[CZTS_host[2]+energy_broadening,CZTS_doping[type][i][7][1]]],
                                                                [[-energy_broadening,CZTS_doping[type][j][7][0]],[CZTS_host[2]+energy_broadening,CZTS_doping[type][j][7][1]]]]
                                                #look for intersection and write in table ["type",[i,i+j+1],x_intersection,y_intersection]
                                                ionisationenergies.append(line_intersection(lines_intersect[0],lines_intersect[1]))
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],14+type-0.41,14+type+0.41,color='black',lw=3.5)
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],14+type-0.4,14+type+0.4,color=color_ionisation(line_intersection(lines_intersect[0],lines_intersect[1])[1]/4),lw=3)
                                                if type < 3:
                                                        ax_CZTS_S.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'*', markerfacecolor='none', markersize=8, markeredgecolor=plot_color[14+type]) #Marker at intersection
                                                else:
                                                        ax_CZTS_I.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'*', markerfacecolor='none', markersize=8, markeredgecolor=plot_color[14+type]) #Marker at intersection
                                else :
                                        if i != j and i != i+j :
                                                ionisationenergies.append([i,i+j])
                                                #look for interesction between two consecutive charged defects
                                                lines_intersect = [[[-energy_broadening,CZTS_doping[type][i][7][0]],[CZTS_host[2]+energy_broadening,CZTS_doping[type][i][7][1]]],
                                                                [[-energy_broadening,CZTS_doping[type][j][7][0]],[CZTS_host[2]+energy_broadening,CZTS_doping[type][j][7][1]]]]
                                                #look for intersection and write in table ["type",[i,i+j+1],x_intersection,y_intersection]
                                                ionisationenergies.append(line_intersection(lines_intersect[0],lines_intersect[1]))
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],14+type-0.41,14+type+0.41,color='black',lw=3.5)
                                                ax_CZTS_ionization.hlines(line_intersection(lines_intersect[0],lines_intersect[1])[0],14+type-0.4,14+type+0.4,color=color_ionisation(line_intersection(lines_intersect[0],lines_intersect[1])[1]/4),lw=3)
                                                if type < 3:
                                                        ax_CZTS_S.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'*', markerfacecolor='none', markersize=8, markeredgecolor=plot_color[14+type]) #Marker at intersection
                                                else :
                                                        ax_CZTS_I.plot(line_intersection(lines_intersect[0],lines_intersect[1])[0],line_intersection(lines_intersect[0],lines_intersect[1])[1],'*', markerfacecolor='none', markersize=8, markeredgecolor=plot_color[14+type]) #Marker at intersection

        ax_CZTS_ionization.hlines(y=CZTS_host[2], xmin=-defect_number, xmax=defect_number, color="darkorange", lw=0.75, label="Conduction band") #conduction band
        ax_CZTS_ionization.fill_between([-defect_number,defect_number],CZTS_host[2],CZTS_host[2]+10,color="darkorange", alpha=0.1) #conduction band
        ax_CZTS_ionization.hlines(y=0, xmin=-defect_number, xmax=defect_number, color="navy", lw=0.75, label="Valence band") #valence band
        ax_CZTS_ionization.fill_between([-defect_number,defect_number],-10,0,color="navy", alpha=0.1) #valence band
        ax_CZTS_ionization.tick_params(axis='y',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=12)
        ax_CZTS_ionization.tick_params(axis='x',length=0,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=11)
        ax_CZTS_ionization.set_ylabel("$\epsilon_{q,q'}(\\alpha)$ [eV]",fontdict = fontlabel , labelpad=5)
        ax_CZTS_ionization.set_xlabel("Defect $\\alpha$",fontdict = fontlabel , labelpad=5)
        ax_CZTS_ionization.set_xlim(-0.5,defect_number-0.5)
        ax_CZTS_ionization.set_xticks(np.arange(0,defect_number,1))
        ax_CZTS_ionization.set_ylim(-0.5,CZTS_host[2]+0.5)
        ax_CZTS_ionization.set_yticks(np.arange(-0.6,CZTS_host[2]+0.6,0.2))
        ax_CZTS_ionization.set_xticklabels(CZTS_Defect_label)
        for i in range(defect_number):
                ax_CZTS_ionization.vlines(x=i-0.5, ymin=-10, ymax=10, color="black", lw=0.25,linestyle='dashed')
        colormap = sns.color_palette("magma", as_cmap=True)
        normalize = mcolors.Normalize(vmin=0, vmax=3)
        s_map = plt.cm.ScalarMappable(norm=normalize, cmap=colormap)
        s_map.set_array([0,1,2,3])
        cbar = fig_CZTS_ionization.colorbar(s_map, ax=ax_CZTS_ionization, spacing='proportional')
        cbarlabel = r'$\beta$ [eV]'
        cbar.set_label(cbarlabel, font=fontlabelbar)

CZTS_transition(mu_cu,mu_zn,mu_sn,CZTS_host[1],Madelung)

fig_CZTS_path = plt.figure(figsize=(8,5))
fig_CZTS_path.subplots_adjust(left=0.129,bottom=0.171,right=0.85,top=0.97)
ax_CZTS_path = plt.subplot()

#-----------------------------------------------------------------------------------------------------------------------------#
#                                                   Equilibrium fermi value                                                   #
#-----------------------------------------------------------------------------------------------------------------------------#

def findFermivalueCZTS(mu_cu,mu_zn,mu_sn,VBM,Madelung):
        Fermienergy_out = 0
        T = 300
        meffe = 0.185
        meffh = 0.96
        EG = 1.5
        Npossible = [["VCu",16],["VSn",8],["VZn",8],["VS",32],
                     ["CuZn",8],["ZnCu",16],["CuSn",8],["SnCu",16],["ZnSn",8],["SnZn",8],
                     ["CuI",100],["SnI",100],["ZnI",100],["Si",100]]
        NC = 2*m.pow((2*m.pi*meffe*cst.electron_mass*cst.Boltzmann*T)/(cst.Planck*cst.Planck),3/2)/1E6
        NV = 2*m.pow((2*m.pi*meffh*cst.electron_mass*cst.Boltzmann*T)/(cst.Planck*cst.Planck),3/2)/1E6
        nEfvalues = 1000
        Concentrationalpha = np.zeros(nEfvalues)
        index = 0
        for EF in np.linspace(0,EG*nEfvalues,nEfvalues)/nEfvalues:
                for defect in range(len(CZTS_Defect)):
                        energy_out = 10
                        chargeindex = 0
                        # Formation energy at EF
                        for charge,chargevalue in enumerate(CZTS_Defect[defect]): # Defect charge state
                                energy_value = CZTS_Defect_calc_all(mu_cu,mu_zn,mu_sn,VBM,EF,defect,charge,Madelung) # Calculate the formation energy
                                if energy_value < energy_out:
                                        energy_out = energy_value
                                        chargeindex = charge
                        Concentrationalpha[index] +=  -CZTS_Defect[defect][chargeindex][1] * Npossible[defect][1] * m.exp(-(energy_out*cst.elementary_charge)/(cst.Boltzmann*T)) #for 64 atoms
                Concentrationalpha[index] = (Concentrationalpha[index]/1309.763788)*1E24 #1cm3 atoms
                index += 1

        index = 0
        for EF in np.linspace(0,EG*nEfvalues,nEfvalues)/nEfvalues:
                n = NC * m.exp(-(EG*cst.elementary_charge - EF*cst.elementary_charge)/(cst.Boltzmann*T))
                p = NV * m.exp(-(EF*cst.elementary_charge)/(cst.Boltzmann*T))
                if abs(abs(p-n) - abs(Concentrationalpha[index])) < 5E14 and np.sign(p-n) == np.sign(Concentrationalpha[index]):
                        #print(EF)
                        Fermienergy_out = EF
                index += 1
        table = np.linspace(0,EG*nEfvalues,nEfvalues)/nEfvalues
        Fermienergy_out = table[np.argmin(abs(Concentrationalpha))]
        return Fermienergy_out

findFermivalueCZTS(mu_cu,mu_zn,mu_sn,CZTS_host[1],Madelung)

#-----------------------------------------------------------------------------------------------------------------------------#
#                                                           Path plot                                                         #
#-----------------------------------------------------------------------------------------------------------------------------#

fig_Fermievolution = plt.figure(figsize=(8,5))
ax_Fermi = plt.subplot()

def CZTS_E_F_path(VBM,Madelung):
        CZTS_chemical_path = np.loadtxt(chemicalpotentialpath_CZTS,delimiter=",") #mu_Cu mu_Zn mu_Sn mu_S
        CZTS_chemical = [0,0,0,0] #mu_Cu mu_Zn mu_Sn mu_S
        CZTS_chemical_save = [0,0,0,0]
        Fermiequilibriumvalue = np.zeros(len(CZTS_chemical_path[0]))

        xpoint = 0
        print("Path plot")
        for point in range(len(CZTS_chemical_path[0])): # for a given point in the phase diagram 1) calculate chemical potentials, 2) Defect formation energy at 0 and EF, 3) ionisation energies

                CZTS_chemical = [CZTS_chemical_path[0][point],CZTS_chemical_path[1][point],CZTS_chemical_path[2][point],0] #mu_Cu mu_Zn mu_Sn mu_S
                CZTS_chemical[3] = (Formation_CZTS - 2*CZTS_chemical_path[0][point] - CZTS_chemical_path[1][point] - CZTS_chemical_path[2][point])/4 #mu_S
                #print('mu_S',CZTS_chemical[3])
                Fermiequilibriumvalue[point] = findFermivalueCZTS(CZTS_chemical[0],CZTS_chemical[1],CZTS_chemical[2],VBM,Madelung) #Compute the equilibrium fermi level for the given concentration.
                #print(Fermiequilibriumvalue[point])

                if CZTS_chemical != CZTS_chemical_save:
                        for defect in range(len(CZTS_Defect)): # Defect type 
                                # Compute formation energy for one defect type, in each charge state and output the minimal formation energy for each fermi energy value, then plot it
                                energy_out = 10
                                chargeindex = 0
                                # Formation energy at EF
                                for charge in range(len(CZTS_Defect[defect])): # Defect charge state
                                        energy_value = CZTS_Defect_calc_all(CZTS_chemical[0],CZTS_chemical[1],CZTS_chemical[2],VBM,Fermiequilibriumvalue[point],defect,charge,Madelung) # Calculate the formation energy
                                        if energy_value < energy_out:
                                                energy_out = energy_value
                                                chargeindex = charge

                                if xpoint == 4 and energy_out < 3:
                                        if defect <= 3:
                                                ax_CZTS_path.plot(xpoint,energy_out,'o', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect],label=CZTS_Defect[defect][chargeindex][0])
                                        elif defect > 3 and defect <=9:
                                                ax_CZTS_path.plot(xpoint,energy_out,'s', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect],label=CZTS_Defect[defect][chargeindex][0])
                                        elif defect > 9:
                                                ax_CZTS_path.plot(xpoint,energy_out,'D', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect],label=CZTS_Defect[defect][chargeindex][0])
                                else:
                                        if defect <= 3:
                                                ax_CZTS_path.plot(xpoint,energy_out,'o', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect])
                                        elif defect > 3 and defect <=9:
                                                ax_CZTS_path.plot(xpoint,energy_out,'s', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect])
                                        elif defect > 9:
                                                ax_CZTS_path.plot(xpoint,energy_out,'D', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect])
                        
                        for defect in range(len(CZTS_doping)): # Defect type 
                                # Compute formation energy for one defect type, in each charge state and output the minimal formation energy for each fermi energy value, then plot it
                                energy_out = 10
                                chargeindex = 0
                                # Formation energy at EF
                                for charge in range(len(CZTS_doping[defect])): # Defect charge state
                                        energy_value = CZTS_Defect_doping(CZTS_chemical[0],CZTS_chemical[1],CZTS_chemical[2],VBM,Fermiequilibriumvalue[point],defect,charge,Madelung) # Calculate the formation energy
                                        if energy_value < energy_out:
                                                energy_out = energy_value
                                                chargeindex = charge

                                if xpoint == 4 and energy_out < 3:
                                        ax_CZTS_path.plot(xpoint,energy_out,'*', markerfacecolor='none', markersize=6, markeredgecolor=plot_color[14+defect],label=CZTS_doping[defect][chargeindex][0])
                                else:
                                        ax_CZTS_path.plot(xpoint,energy_out,'*', markerfacecolor='none', markersize=6, markeredgecolor=plot_color[14+defect])
                        if xpoint == 4:
                                ax_Fermi.plot(xpoint,Fermiequilibriumvalue[point],'.', markerfacecolor='none', markersize=6, markeredgecolor='#006591',label="$Cu_2ZnSnS_4$")
                        else:
                                ax_Fermi.plot(xpoint,Fermiequilibriumvalue[point],'.', markerfacecolor='none', markersize=6, markeredgecolor='#006591')
                        xpoint += 1

                CZTS_chemical_save = CZTS_chemical

        ax_CZTS_path.hlines(y=0, xmin=-10, xmax=100, color="black", lw=0.25,linestyle='dashed')                        # 0 formation energy line
        ax_CZTS_path.set_ylabel("$\Delta H_F$ [eV]",fontdict = fontlabel , labelpad=5)
        ax_CZTS_path.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)
        ax_CZTS_path.legend(loc='upper left',prop={'size': 12},frameon=False,bbox_to_anchor=(1.0, 1))
        ax_CZTS_path.set_xlabel("$\mu_i$ [eV]",fontdict = fontlabel , labelpad=5)
        xtick=np.linspace(0,len(CZTS_chemical_path[0])-8,9)
        ax_CZTS_path.set_xticks(xtick)
        ax_CZTS_path.set_xticklabels(["A", "B", "C", "D", "E", "F","G","H","I"])
        ax_CZTS_path.set_xlim(-1,len(CZTS_chemical_path[0])-8+1)
        ax_CZTS_path.set_ylim(-0.5,3)
        ax_CZTS_path.set_yticks(np.arange(-1,7,1)/2)

        ax_Fermi.set_ylabel("$E_{F,eq}$ [eV]",fontdict = fontlabel , labelpad=5)
        ax_Fermi.set_xticks(xtick)
        ax_Fermi.set_xticklabels(["A", "B", "C", "D", "E", "F","G","H","I"])
        ax_Fermi.legend(prop={'size': 12})
        ax_Fermi.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)

CZTS_E_F_path(CZTS_host[1],Madelung)

#-----------------------------------------------------------------------------------------------------------------------------#
#                                                       Sliders                                                               #
#-----------------------------------------------------------------------------------------------------------------------------#

if Set_Slider:
        axCu = fig_CZTS_all.add_axes([0.1, 0.17, 0.35, 0.02], facecolor='orange')
        cu_slide = Slider(axCu,"$\mu_{Cu}$", -2, 0, valinit=mu_cu, valstep=0.001,color='lightgrey')
        cu_slide.label.set_size(14)
        axZn = fig_CZTS_all.add_axes([0.1, 0.12, 0.35, 0.02], facecolor='orange')
        zn_slide = Slider(axZn,"$\mu_{Zn}$", -2, 0, valinit=mu_zn, valstep=0.001,color='lightgrey')
        zn_slide.label.set_size(14)
        axSn = fig_CZTS_all.add_axes([0.1, 0.07, 0.35, 0.02], facecolor='orange')
        sn_slide = Slider(axSn,"$\mu_{Sn}$", -2, 0, valinit=mu_sn, valstep=0.001,color='lightgrey')
        sn_slide.label.set_size(14)

        def Update_Defect_Plot_CZTS(val):
                mu_cu = cu_slide.val
                mu_zn = zn_slide.val
                mu_sn = sn_slide.val
                fermi = CZTS_host[1]
                ax_CZTS_V.cla()
                ax_CZTS_S.cla()
                ax_CZTS_I.cla()
                CZTS_Defect_plot(mu_cu,mu_zn,mu_sn,fermi,Madelung)
                fig_CZTS_all.canvas.draw_idle()

        cu_slide.on_changed(Update_Defect_Plot_CZTS)
        zn_slide.on_changed(Update_Defect_Plot_CZTS)
        sn_slide.on_changed(Update_Defect_Plot_CZTS)

if savefigs == True:
        fig_CZTS_all.savefig(pathsave + 'CZTS_EF.pdf')
        fig_CZTS_ionization.savefig(pathsave + 'CZTS_ionization.pdf')
        fig_CZTS_path.savefig(pathsave + 'CZTS_path.pdf')

######################################################
#####                                           ###### 
#####                    CZGS                   ###### 
#####                                           ###### 
######################################################

if CZGS_calc:
        ########################################################################################################################
        #                                                                                                                      #
        #                                           Formation energy plot for CZGS                                             #
        #                                                                                                                      #
        ########################################################################################################################

        Pure_energies_CZGS = [["Cu",-14.58743314,-3.65226],["Zn",-2.53611365,-1.2574],["Ge",-36.27147174,-5.35416],["S",-167.63710169,-5.23831]]
        Formation_CZGS = 8*(-315.53754352/64) - 2 * Pure_energies_CZGS[0][2] - Pure_energies_CZGS[1][2] - Pure_energies_CZGS[2][2] - 4*Pure_energies_CZGS[3][2]

        #[Energy, VBM, Eg]
        CZGS_host = [-315.53754352,3.84721361,1.891] #Fermi in DOSCAR CZGS supercell = 3.84721361 eV Eg = 1.720618
        #Fermi in DOSCAR CZGS 8 atoms = 3.611929 eV Eg = 1.720618

        #[0-"Name", 1-charge state, 2-n_Cu, 3-n_Zn, 4-n_X, 5-n_S, 6-total energy, 7-Formation energy, 8-Potential correction compare to host supercell]                     
        CZGS_Defect = [[["$V_{Cu}$",0,-1,0,0,0,-310.367005745,[0,0],0.09729337],["$V_{Cu}^{-1}$",-1,-1,0,0,0,-306.946003509,[0,0],0.17214677]],
                [["$V_{Ge}$",0,0,0,-1,0,-304.368076735,[0,0],0.17402669],["$V_{Ge}^{-1}$",-1,0,0,-1,0,-300.766257792,[0,0],0.27967037],["$V_{Ge}^{-2}$",-2,0,0,-1,0,-297.230356411,[0,0],0.38444063],["$V_{Ge}^{-3}$",-3,0,0,-1,0,-293.71363683,[0,0],0.48503162],["$V_{Ge}^{-4}$",-4,0,0,-1,0,-290.239280589,[0,0],0.5783622]],
                [["$V_{Zn}$",0,0,-1,0,0,-310.62843,[0,0],0.04340232],["$V_{Zn}^{-1}$",-1,0,-1,0,0,-307.079941104,[0,0],0.12488838],["$V_{Zn}^{-2}$",-2,0,-1,0,0,-303.679692105,[0,0],0.20709505]],
                [["$V_S$",0,0,0,0,-1,-308.064242701,[0,0],0.18994073],["$V_S^{+1}$",+1,0,0,0,-1,-311.712705259,[0,0],0.09894331]],
                [["$Cu_{Zn}$",0,+1,-1,0,0,-316.266403777,[0,0],-0.0471987],["$Cu_{Zn}^{-1}$",-1,+1,-1,0,0,-312.507872086,[0,0],0.02372229]],
                [["$Zn_{Cu}$",0,-1,+1,0,0,-312.410550313,[0,0],0.0202166],["$Zn_{Cu}^{+1}$",+1,-1,+1,0,0,-318.280393211,[0,0],-0.19022979]],
                [["$Cu_{Ge}$",0,+1,0,-1,0,-310.866792681,[0,0],0.06295009],["$Cu_{Ge}^{-1}$",-1,+1,0,-1,0,-307.157044774,[0,0],0.14331469],["$Cu_{Ge}^{-2}$",-2,+1,0,-1,0,-303.297852351,[0,0],0.23275869],["$Cu_{Ge}^{-3}$",-3,+1,0,-1,0,-299.448017363,[0,0],0.33638823]],
                [["$Ge_{Cu}$",0,-1,0,+1,0,-314.041462885,[0,0],-0.00970278],["$Ge_{Cu}^{+1}$",+1,-1,0,+1,0,-319.871295789,[0,0],-0.0135146],["$Ge_{Cu}^{+2}$",+2,-1,0,+1,0,-324.043767392,[0,0],-0.00239415],["$Ge_{Cu}^{+3}$",+3,-1,0,+1,0,-329.085757356,[0,0],0.25928952]],
                [["$Zn_{Ge}$",0,0,+1,-1,0,-310.243531271,[0,0],0.11996433],["$Zn_{Ge}^{-1}$",-1,0,+1,-1,0,-306.719258439,[0,0],0.19587446],["$Zn_{Ge}^{-2}$",-2,0,+1,-1,0,-303.279328688,[0,0],0.27003811]],
                [["$Ge_{Zn}$",0,0,-1,+1,0,-316.935506119,[0,0],-0.01704416],["$Ge_{Zn}^{+1}$",+1,0,-1,+1,0,-321.54176634,[0,0],-0.04440281],["$Ge_{Zn}^{+2}$",+2,0,-1,+1,0,-326.985282072,[0,0],0.17321663]],
                [["$Cu_{I}$",0,+1,0,0,0,-316.925615018,[0,0],-0.04834668],["$I_{Cu}^{+1}$ ",+1,+1,0,0,0,-322.791110519,[0,0],-0.08583074]],
                [["$Ge_{I}$",0,0,0,+1,0,-316.966350453,[0,0],-0.03057512],["$I_{Ge}^{+1}$",+1,0,0,+1,0,-321.311164798,[0,0],-0.09535345],["$I_{Ge}^{+2}$",+2,0,0,+1,0,-327.008605473,[0,0],-0.01892878],["$I_{Ge}^{+3}$",+3,0,0,+1,0,-331.08329483,[0,0],-0.07429894],["$I_{Ge}^{+4}$",+4,0,0,+1,0,-335.274470985,[0,0],-0.0885865]],
                [["$Zn_{I}$",0,0,+1,0,0,-313.620057874,[0,0],0.06640339],["$I_{Zn}^{+1}$",+1,0,+1,0,0,-318.769631174,[0,0],0.07516314],["$I_{Zn}^{+2}$",+2,0,+1,0,0,-324.761627955,[0,0],-0.0339012]],
                [["$S_I$",0,0,0,0,+1,-316.857706301,[0,0],0.37170205],["$I_S^{+1}$",-1,0,0,0,+1,-311.830049233,[0,0],0.10397123]]
                ]

        #C_sh, Epsilon, Alpha_M, L
        Madelung_CZGS = [-0.342,6.44*cst.epsilon_0,Madelung,m.pow(1225.45*m.pow(10, -30),1/3)] #1.63664

        fig_CZGS_all, (ax_CZGS_V, ax_CZGS_S, ax_CZGS_I) = plt.subplots(1,3,figsize=(12,7))
        fig_CZGS_all.subplots_adjust(bottom=0.3, wspace = 0.2,right=0.98,left=0.1,top=0.98)
        plot_color = ["#00708F","#E6A100","#E61F00","#4D009B","#B98000","#00A053","#00BEB6","#005228","#754D00","#720008","#FF4365","#2B50AA","#FF9FE5","#73A6AD"]
        subplotarry_CZGS = [ax_CZGS_V,ax_CZGS_S,ax_CZGS_I]

        #-----------------------------------------------------------------------------------------------------------------------------#
        #                                               Formation Energy Calculations                                                 #
        #-----------------------------------------------------------------------------------------------------------------------------#

        def CZGS_Defect_calc_all(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,Fermi,x,type,charge,Madelung): # Function to compute the defect formation energy at a given energy "x" (Fermi position)
                CZGS_chemical = [CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,0] #mu_Cu mu_Zn mu_Ge
                CZGS_chemical[3] = (Formation_CZGS - 2*CZGS_chemical[0] - CZGS_chemical[1] - CZGS_chemical[2])/4 #mu_S
                out = 0
                #Hf    #E(i,h)                       #Ehost         #charge                       #VBM         #Potentiel alignement
                out = CZGS_Defect[type][charge][6] - CZGS_host[0] + CZGS_Defect[type][charge][1]*(Fermi + x) - CZGS_Defect[type][charge][1]*CZGS_Defect[type][charge][8] + (0.65)*((m.pow(CZGS_Defect[type][charge][1]*cst.elementary_charge,2)*Madelung)/(2*Madelung_CZGS[1]*Madelung_CZGS[3]))/cst.elementary_charge
                for j in range(4):
                        #Hf          #ni                             #mu_i             #E_i
                        out = out - CZGS_Defect[type][charge][2+j]*(CZGS_chemical[j] + Pure_energies_CZGS[j][2])
                return out

        #-----------------------------------------------------------------------------------------------------------------------------#
        #                                              Defect Plot Function                                                           #
        #-----------------------------------------------------------------------------------------------------------------------------#

        def CZGS_Defect_plot(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,Madelung): # Plot function for defect formation energy V3
                init = -1 # EF init value (within valence band)
                end = CZGS_host[2] + 1 # EF end value (within conduction band)
                x = np.linspace(init,end,nline) #1000 pts from init to end for each line
                
                print("Defect plot")
                for defect,defectvalue in enumerate(CZGS_Defect): # Defect type 
                        # Compute formation energy for one defect type, in each charge state and output the minimal formation energy for each fermi energy value, then plot it
                        line_defect_out = np.zeros(len(x))
                        line_defect = np.zeros(len(x))
                        line_defect_temp = np.zeros(len(x))

                        # Formation energy at a given energy E
                        for charge,chargevalue in enumerate(CZGS_Defect[defect]): # Defect charge state
                                for a in range(len(x)): # Loop over the energies
                                        line_defect_temp[a] = CZGS_Defect_calc_all(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,x[a],defect,charge,Madelung) # Calculate the formation energy
                                line_defect = np.vstack((line_defect,line_defect_temp)) # Put 2 lines arrays in a same arrays as [[...,...],[...,...]]
                        line_defect = np.delete(line_defect, 0, 0) # Remove first line which is all 0 elements 
                        # Find the minimum of all lines (i.e. minimum of formation energy for each charge state of one defect) and also gives the energy at which the deviation occurs
                        line_defect_out = np.min(line_defect, axis=0)

                        # Plot of the line
                        if(defect <= 3):    
                                ax_CZGS_V.plot(x,line_defect_out,lw = 1, label=CZGS_Defect[defect][0][0], color=plot_color[defect])
                                for charge in range(len(CZGS_Defect[defect])): # Defect charge state
                                        if draweachdefect == True:
                                                ax_CZGS_V.plot(x,line_defect[charge],lw = 0.5, color=plot_color[defect],linestyle='dashed')
                        elif(defect > 3 and defect <=9):
                                ax_CZGS_S.plot(x,line_defect_out,lw = 1, label=CZGS_Defect[defect][0][0], color=plot_color[defect])
                                for charge in range(len(CZGS_Defect[defect])): # Defect charge state
                                        if draweachdefect == True:
                                                ax_CZGS_S.plot(x,line_defect[charge],lw = 0.5, color=plot_color[defect],linestyle='dashed')
                        elif(defect > 9):
                                ax_CZGS_I.plot(x,line_defect_out,lw = 1, label=CZGS_Defect[defect][0][0], color=plot_color[defect])
                                for charge in range(len(CZGS_Defect[defect])): # Defect charge state
                                        if draweachdefect == True:
                                                ax_CZGS_I.plot(x,line_defect[charge],lw = 0.5, color=plot_color[defect],linestyle='dashed')

                for defectplot in subplotarry_CZGS:
                        defectplot.vlines(x=0, ymin=-10, ymax=10, color="navy", lw=0.75)                                            # label="$\epsilon_{VBM}$"
                        defectplot.vlines(x=CZGS_host[2], ymin=-10, ymax=10, color="darkorange", lw=0.75)                           # label="$E_{G}$"
                        defectplot.vlines(x=(nbkb_CZGS*cst.Boltzmann*300)/cst.elementary_charge, ymin=-10, ymax=10, color="maroon",lw=1.5, linestyle = '--') 
                        defectplot.hlines(y=0, xmin=-10, xmax=10, color="black", lw=0.25,linestyle='dashed')                        # 0 formation energy line
                        defectplot.fill_between(np.arange(-10,0.01,0.01),-10,10,color="navy", alpha=0.1)                            # Conduction band
                        defectplot.fill_between(np.arange(CZGS_host[2],CZGS_host[2]+5,0.01),-10,10,color="darkorange", alpha=0.1)   # Valence band

                ax_CZGS_V.set_ylabel("$\Delta H_F$ [eV]",fontdict = fontlabel , labelpad=5)
                for defectplot in subplotarry_CZGS:
                        defectplot.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)
                        defectplot.legend(loc='upper right',prop={'size': 14})
                        defectplot.set_xlabel("$E_F$ [eV]",fontdict = fontlabel , labelpad=5)
                        defectplot.set_xlim(-0.5,CZGS_host[2]+0.5)
                        defectplot.set_xticks(np.arange(-0.5,2.5,0.5))
                        defectplot.set_ylim(-0.5,yenergy)
                        defectplot.set_yticks(np.arange(-0.5,yenergy+0.5,0.5))

        CZGS_Defect_plot(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,CZGS_host[1],Madelung)  

        #-----------------------------------------------------------------------------------------------------------------------------#
        #                                                 Ionization  Levels                                                          #
        #-----------------------------------------------------------------------------------------------------------------------------#
        #The ionization levels is defined as the the Fermi energy for which a defect formation energy of a defect "a->CZGS_Defect[i]"
        #in a charge state "q->CZGS_Defect[i][j]" is equal to the defect formation energy of the same defect "a->CZGS_Defect[i]" in 
        #another charge state q' -> use intersection function 
        defect_number_CZGS = len(CZGS_Defect)
        energy_broadening = 0 

        fig_CZGS_ionization = plt.figure(figsize=(8,5))
        fig_CZGS_ionization.subplots_adjust(left=0.12,bottom=0.171,right=1,top=0.97)
        ax_CZGS_ionization = plt.subplot()
        CZGS_Defect_label = ["" for x in range(len(CZGS_Defect))] # Set defect name
        for i in range(len(CZGS_Defect)):
                CZGS_Defect_label[i] = CZGS_Defect[i][0][0]

        ionisationenergies_CZGS = []
        def CZGS_transition(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,Madelung):
                print("Ionisation plot")
                color_ionisation = sns.color_palette("magma", as_cmap=True)
                ionisationenergies_CZGS.clear()
                CZGS_chemical = [CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,0] #mu_Cu mu_Zn mu_Sn mu_S
                CZGS_chemical[3] = (Formation_CZGS - 2*CZGS_chemical[0] - CZGS_chemical[1] - CZGS_chemical[2])/4 #mu_S
                for type,typevalue in enumerate(CZGS_Defect): # Defect type
                        for charge,chargevalue in enumerate(CZGS_Defect[type]): # Defect charge state
                                CZGS_Defect[type][charge][7][0] = CZGS_Defect_calc_all(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,0,type,charge,Madelung)
                                CZGS_Defect[type][charge][7][1] = CZGS_Defect_calc_all(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,CZGS_host[2],type,charge,Madelung)
                for type,typevalue in enumerate(CZGS_Defect):
                        ionisationenergies_CZGS.append(CZGS_Defect[type][0][0])
                        #Search for intersection between defects formation energies:
                        for i,ivalue in enumerate(CZGS_Defect[type]):
                                for j,jvalue in enumerate(CZGS_Defect[type]):
                                        if oneelectrontransitiononly == True:
                                                if i != j and i != i+j and abs(i-j) == 1: #last argument is to take into account only 1e transition
                                                        ionisationenergies_CZGS.append([i,i+j])
                                                        #look for interesction between two consecutive charged defects
                                                        lines_intersect_CZGS = [[[-energy_broadening,CZGS_Defect[type][i][7][0]],[CZGS_host[2]+energy_broadening,CZGS_Defect[type][i][7][1]]],
                                                                        [[-energy_broadening,CZGS_Defect[type][j][7][0]],[CZGS_host[2]+energy_broadening,CZGS_Defect[type][j][7][1]]]]
                                                        #look for intersection and write in table ["type",[i,i+j+1],x_intersection,y_intersection]
                                                        ionisationenergies.append(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1]))
                                                        ax_CZGS_ionization.hlines(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],type-0.41,type+0.41,color='black',lw=3.5)
                                                        ax_CZGS_ionization.hlines(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],type-0.4,type+0.4,color=color_ionisation(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1]/4),lw=3)
                                                        if type <= 3:
                                                                ax_CZGS_V.plot(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1],'o', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                        elif type > 3 and type <= 9:
                                                                ax_CZGS_S.plot(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1],'s', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                        elif type > 9:
                                                                ax_CZGS_I.plot(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1],'D', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                        else:
                                                if i != j and i != i+j :
                                                        ionisationenergies_CZGS.append([i,i+j])
                                                        #look for interesction between two consecutive charged defects
                                                        lines_intersect_CZGS = [[[-energy_broadening,CZGS_Defect[type][i][7][0]],[CZGS_host[2]+energy_broadening,CZGS_Defect[type][i][7][1]]],
                                                                        [[-energy_broadening,CZGS_Defect[type][j][7][0]],[CZGS_host[2]+energy_broadening,CZGS_Defect[type][j][7][1]]]]
                                                        #look for intersection and write in table ["type",[i,i+j+1],x_intersection,y_intersection]
                                                        ionisationenergies.append(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1]))
                                                        ax_CZGS_ionization.hlines(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],type-0.41,type+0.41,color='black',lw=3.5)
                                                        ax_CZGS_ionization.hlines(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],type-0.4,type+0.4,color=color_ionisation(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1]/4),lw=3)
                                                        if type <= 3:
                                                                ax_CZGS_V.plot(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1],'o', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                        elif type > 3 and type <= 9:
                                                                ax_CZGS_S.plot(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1],'s', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection
                                                        elif type > 9:
                                                                ax_CZGS_I.plot(line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[0],line_intersection(lines_intersect_CZGS[0],lines_intersect_CZGS[1])[1],'D', markerfacecolor='none', markersize=7, markeredgecolor=plot_color[type]) #Marker at intersection

                ax_CZGS_ionization.hlines(y=CZGS_host[2], xmin=-defect_number_CZGS, xmax=defect_number_CZGS, color="darkorange", lw=0.75, label="Conduction band") #conduction band
                ax_CZGS_ionization.fill_between([-defect_number_CZGS,defect_number_CZGS],CZGS_host[2],CZGS_host[2]+10,color="darkorange", alpha=0.1) #conduction band
                ax_CZGS_ionization.hlines(y=0, xmin=-defect_number_CZGS, xmax=defect_number_CZGS, color="navy", lw=0.75, label="Valence band") #valence band
                ax_CZGS_ionization.fill_between([-defect_number_CZGS,defect_number_CZGS],-10,0,color="navy", alpha=0.1) #valence band
                ax_CZGS_ionization.tick_params(axis='y',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=12)
                ax_CZGS_ionization.tick_params(axis='x',length=0,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=12)
                ax_CZGS_ionization.set_ylabel("$\epsilon_{q,q'}(\\alpha)$ [eV]",fontdict = fontlabel , labelpad=5)
                ax_CZGS_ionization.set_xlabel("Defect $\\alpha$",fontdict = fontlabel , labelpad=5)
                ax_CZGS_ionization.set_xlim(-0.5,defect_number_CZGS-0.5)
                ax_CZGS_ionization.set_xticks(np.arange(0,defect_number_CZGS,1))
                ax_CZGS_ionization.set_ylim(-0.5,CZGS_host[2]+0.5)
                ax_CZGS_ionization.set_yticks(np.arange(-0.6,CZGS_host[2]+0.6,0.2))
                ax_CZGS_ionization.set_xticklabels(CZGS_Defect_label)
                for i in range(defect_number_CZGS):
                        ax_CZGS_ionization.vlines(x=i-0.5, ymin=-10, ymax=10, color="black", lw=0.25,linestyle='dashed')

                colormap = sns.color_palette("magma", as_cmap=True)
                normalize = mcolors.Normalize(vmin=0, vmax=3)
                s_map = plt.cm.ScalarMappable(norm=normalize, cmap=colormap)
                s_map.set_array([0,1,2,3])
                cbar = fig_CZGS_ionization.colorbar(s_map, ax=ax_CZGS_ionization, spacing='proportional')
                cbarlabel = r'$\beta$ [eV]'
                cbar.set_label(cbarlabel, font=fontlabelbar)

        CZGS_transition(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,CZGS_host[1],Madelung)

        fig_CZGS_path = plt.figure(figsize=(8,5))
        fig_CZGS_path.subplots_adjust(left=0.129,bottom=0.171,right=0.85,top=0.97)
        ax_CZGS_path = plt.subplot()

        #-----------------------------------------------------------------------------------------------------------------------------#
        #                                                Fermi equilibrium function                                                   #
        #-----------------------------------------------------------------------------------------------------------------------------#

        def findFermivalueCZGS(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,Madelung):
                Fermienergy_out_CZGS = 0
                T = 300
                meffe = 0.225
                meffh = 1
                EG = 1.85
                Npossible = [["VCu",16],["VGe",8],["VZn",8],["VS",32],
                        ["CuZn",8],["ZnCu",16],["CuGe",8],["GeCu",16],["ZnGe",8],["GeZn",8],
                        ["CuI",100],["GeI",100],["ZnI",100],["Si",100]]
                NC = 2*m.pow((2*m.pi*meffe*cst.electron_mass*cst.Boltzmann*T)/(cst.Planck*cst.Planck),3/2)/1E6
                NV = 2*m.pow((2*m.pi*meffh*cst.electron_mass*cst.Boltzmann*T)/(cst.Planck*cst.Planck),3/2)/1E6
                #print(NC," ", NV)
                nEfvalues = 1000
                Concentrationalpha = np.zeros(nEfvalues)
                index = 0
                for EF in np.linspace(0,EG*nEfvalues,nEfvalues)/nEfvalues:
                        for defect in range(len(CZGS_Defect)):
                                energy_out = 10
                                chargeindex = 0
                                # Formation energy at EF
                                for charge in range(len(CZGS_Defect[defect])): # Defect charge state
                                        energy_value = CZGS_Defect_calc_all(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,VBM,EF,defect,charge,Madelung) # Calculate the formation energy
                                        if energy_value < energy_out:
                                                energy_out = energy_value
                                                chargeindex = charge
                                Concentrationalpha[index] +=  -CZGS_Defect[defect][chargeindex][1] * Npossible[defect][1] * m.exp(-(energy_out*cst.elementary_charge)/(cst.Boltzmann*T)) #for 64 atoms
                        Concentrationalpha[index] = (Concentrationalpha[index]/1225.452891)*1E24 #1cm3 atoms
                        index += 1

                index = 0
                for EF in np.linspace(0,EG*nEfvalues,nEfvalues)/nEfvalues:
                        n = NC * m.exp(-(EG*cst.elementary_charge - EF*cst.elementary_charge)/(cst.Boltzmann*T))
                        p = NV * m.exp(-(EF*cst.elementary_charge)/(cst.Boltzmann*T))
                        #print(EF," ",f"{Concentrationalpha[index]:.2e}"," ",f"{p-n:.2e}")
                        if abs(abs(p-n) - abs(Concentrationalpha[index])) < 1E12 and np.sign(p-n) == np.sign(Concentrationalpha[index]):
                                #print(EF)
                                Fermienergy_out_CZGS = EF
                        index += 1
                
                table = np.linspace(0,EG*nEfvalues,nEfvalues)/nEfvalues
                Fermienergy_out_CZGS = table[np.argmin(abs(Concentrationalpha))]

                return Fermienergy_out_CZGS

        findFermivalueCZGS(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,CZGS_host[1],Madelung)

        #-----------------------------------------------------------------------------------------------------------------------------#
        #                                                           Path plot                                                         #
        #-----------------------------------------------------------------------------------------------------------------------------#        

        def CZGS_E_F_path(VBM,Madelung):
                print("Path plot")
                CZGS_chemical_path = np.loadtxt(chemicalpotentialpath_CZGS,delimiter=",") #mu_Cu mu_Zn mu_Ge mu_S
                CZGS_chemical = [0,0,0,0] #mu_Cu mu_Zn mu_Sn mu_S
                CZGS_chemical_save = [0,0,0,0]
                xpoint = 0
                Fermiequilibriumvalue_CZGS = np.zeros(len(CZGS_chemical_path[0]))

                for point in range(len(CZGS_chemical_path[0])): # for a given point in the phase diagram 1) calculate chemical potentials, 2) Defect formation energy at 0 and EF, 3) ionisation energies
                        
                        CZGS_chemical = [CZGS_chemical_path[0][point],CZGS_chemical_path[1][point],CZGS_chemical_path[2][point],0] #mu_Cu mu_Zn mu_Sn mu_S
                        CZGS_chemical[3] = (Formation_CZTS - 2*CZGS_chemical_path[0][point] - CZGS_chemical_path[1][point] - CZGS_chemical_path[2][point])/4 #mu_S
                        #print('mu_s',CZGS_chemical[3])
                        Fermiequilibriumvalue_CZGS[point] = findFermivalueCZGS(CZGS_chemical[0],CZGS_chemical[1],CZGS_chemical[2],VBM,Madelung) #Compute the equilibrium fermi level for the given concentration.

                        if CZGS_chemical != CZGS_chemical_save:

                                for defect in range(len(CZGS_Defect)): # Defect type 
                                        # Compute formation energy for one defect type, in each charge state and output the minimal formation energy for each fermi energy value, then plot it
                                        energy_out = 10
                                        defectindex = 0

                                        # Formation energy at EF
                                        for charge in range(len(CZGS_Defect[defect])): # Defect charge state
                                                energy_value = CZGS_Defect_calc_all(CZGS_chemical[0],CZGS_chemical[1],CZGS_chemical[2],VBM,Fermiequilibriumvalue_CZGS[point],defect,charge,Madelung) # Calculate the formation energy
                                                if energy_value < energy_out:
                                                    energy_out = energy_value
                                                    defectindex = charge
                
                                        if xpoint == 4 and energy_out < 3:
                                                if defect <= 3:
                                                        ax_CZGS_path.plot(xpoint,energy_out,'o', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect],label=CZGS_Defect[defect][defectindex][0])
                                                elif defect > 3 and defect <=9:
                                                        ax_CZGS_path.plot(xpoint,energy_out,'s', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect],label=CZGS_Defect[defect][defectindex][0])
                                                elif defect > 9:
                                                        ax_CZGS_path.plot(xpoint,energy_out,'D', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect],label=CZGS_Defect[defect][defectindex][0])
                                        else:
                                                if defect <= 3:
                                                        ax_CZGS_path.plot(xpoint,energy_out,'o', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect])
                                                elif defect > 3 and defect <=9:
                                                        ax_CZGS_path.plot(xpoint,energy_out,'s', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect])
                                                elif defect > 9:
                                                        ax_CZGS_path.plot(xpoint,energy_out,'D', markerfacecolor='none', markersize=4, markeredgecolor=plot_color[defect])
                                
                                if xpoint == 4:
                                        ax_Fermi.plot(xpoint,Fermiequilibriumvalue_CZGS[point],'.', markerfacecolor='none', markersize=6, markeredgecolor='#E60B00',label="$Cu_2ZnGeS_4$")
                                else:
                                        ax_Fermi.plot(xpoint,Fermiequilibriumvalue_CZGS[point],'.', markerfacecolor='none', markersize=6, markeredgecolor='#E60B00')
                                xpoint += 1
                        CZGS_chemical_save = CZGS_chemical

                ax_CZGS_path.hlines(y=0, xmin=-10, xmax=100, color="black", lw=0.25,linestyle='dashed')                        # 0 formation energy line
                ax_CZGS_path.set_ylabel("$\Delta H_F$ [eV]",fontdict = fontlabel , labelpad=5)
                ax_CZGS_path.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)
                ax_CZGS_path.legend(loc='upper left',prop={'size': 12},frameon=False,bbox_to_anchor=(1.0, 0.9))
                ax_CZGS_path.set_xlabel("$\mu_i$ [eV]",fontdict = fontlabel , labelpad=5)
                xtick=np.linspace(0,len(CZGS_chemical_path[0])-8,9)
                ax_CZGS_path.set_xticks(xtick)
                ax_CZGS_path.set_xticklabels(["A", "B", "C", "D", "E", "F","G","H","I"])
                ax_CZGS_path.set_xlim(-1,len(CZGS_chemical_path[0])-8+1)
                ax_CZGS_path.set_ylim(-0.5,3)
                ax_CZGS_path.set_yticks(np.arange(-1,7,1)/2)
                ax_Fermi.legend(prop={'size': 12})

        CZGS_E_F_path(CZGS_host[1],Madelung)

        #-----------------------------------------------------------------------------------------------------------------------------#
        #                                                       Sliders                                                               #
        #-----------------------------------------------------------------------------------------------------------------------------#

        if Set_Slider:
                CZGS_axCu = fig_CZGS_all.add_axes([0.57, 0.17, 0.35, 0.02], facecolor='orange')
                CZGS_cu_slide = Slider(CZGS_axCu,"$\mu_{Cu}$", -2, 0, valinit=CZGS_mu_cu, valstep=0.001,color='lightgrey')
                CZGS_cu_slide.label.set_size(14)
                CZGS_axZn = fig_CZGS_all.add_axes([0.57, 0.12, 0.35, 0.02], facecolor='orange')
                CZGS_zn_slide = Slider(CZGS_axZn,"$\mu_{Zn}$", -2, 0, valinit=CZGS_mu_zn, valstep=0.001,color='lightgrey')
                CZGS_zn_slide.label.set_size(14)
                CZGS_axGe = fig_CZGS_all.add_axes([0.57, 0.07, 0.35, 0.02], facecolor='orange')
                CZGS_ge_slide = Slider(CZGS_axGe,"$\mu_{Ge}$", -2, 0, valinit=CZGS_mu_ge, valstep=0.001,color='lightgrey')
                CZGS_ge_slide.label.set_size(14)

                def Update_Defect_Plot_CZGS(val):
                        CZGS_mu_cu = CZGS_cu_slide.val
                        CZGS_mu_zn = CZGS_zn_slide.val
                        CZGS_mu_ge = CZGS_ge_slide.val
                        CZGS_VBM = CZGS_host[1]
                        CZGS_madelung = Madelung
                        ax_CZGS_V.cla()
                        ax_CZGS_S.cla()
                        ax_CZGS_I.cla()
                        CZGS_Defect_plot(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,CZGS_VBM,Madelung)
                        fig_CZGS_all.canvas.draw_idle()
                        ax_CZGS_ionization.cla()
                        CZGS_transition(CZGS_mu_cu,CZGS_mu_zn,CZGS_mu_ge,CZGS_VBM,CZGS_madelung)
                        fig_CZGS_ionization.canvas.draw_idle()

                CZGS_cu_slide.on_changed(Update_Defect_Plot_CZGS)
                CZGS_zn_slide.on_changed(Update_Defect_Plot_CZGS)
                CZGS_ge_slide.on_changed(Update_Defect_Plot_CZGS)

        if savefigs == True:
                fig_CZGS_all.savefig(pathsave + 'CZGS_EF.pdf')
                fig_CZGS_ionization.savefig(pathsave + 'CZGS_ionization.pdf')
                fig_CZGS_path.savefig(pathsave + 'CZGS_path.pdf')
                fig_Fermievolution.savefig(pathsave + 'Fermi_path.pdf')

plt.show()
