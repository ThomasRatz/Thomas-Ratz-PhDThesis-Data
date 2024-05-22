                                        #################################################
                                        #                                               #
                                        #       2D plot of Kesterite Phase Diagram      #
                                        #                                               #
                                        #################################################
#Librairies using python 3.9.6

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from celluloid import Camera # getting the camera
#plt.xkcd() #If you want some fun uncomment this and rerun !

#Fonts
#-----
font = {'family': 'serif', 'size': 20}
fontlabel = {'family': 'serif', 'size': 20}
axcolor = 'lightgrey'

#Parameters
#---------
mu_cu_init = -0.273 #[-0.273,-0.55,-0.82] from Cu rich to Cu poor as presented in Fig.5.1
nbr = 1000          #Nbr of points in lines
Gecompound = True  #Compute CZGS phase  diagram
saveflag = False    #Save figures
gif=True            #Generate GIF
path = "/Users/thomasratz/Desktop/"   #"Path to save figures"

#HSE06 input data from DFT calculations
#--------------------------------------

#Pure phases with total energy and total energy per atom [eV/atom] relaxed using SCAN functionnal and total energy obtained by HSE06 electronic relaxation.
Pure_energies_CZTS = [["Cu",-3.65226*4,-3.65226],["Zn",-1.2574*2,-1.2574],["S",-5.23831*32,-5.23831],["Sn",-4.53397*8,-4.53397]]
Pure_energies_CZGS = [["Cu",-3.65226*4,-3.65226],["Zn",-1.2574*2,-1.2574],["S",-5.23831*32,-5.23831],["Ge",-5.35416*8,-5.35416]]

#Secondary phases with total energy and total energy per atom [eV/atom] relaxed using SCAN functionnal and total energy obtained by HSE06 electronic relaxation.
Secondaries_energies_CZTS = [["Cu7S4",-197.83995569,-197.83995569/44],["CuS",-56.30473445,-56.30473445/12],["CuS2",-28.68865111,-28.68865111/6],
                            ["ZnS",-33.56936371,-33.56936371/8],["ZnS2",-52.62629248,-52.62629248/12],["Zn8Cu5",-119.03567336,-119.03567336/52],
                            ["Zn35Cu17",-112.33610199,-112.33610199/52],["Cu2SnS3",-120.19605530,-120.19605530/24],["CuSn",-16.41968107,-16.41968107/4],
                            ["SnS",-42.47698170,-42.47698170/8],["SnS2",-16.24787812,-16.24787812/3]]
Secondaries_energies_CZGS = [["Cu7S4",-197.83995569,-197.83995569/44],["CuS",-56.30473445,-56.30473445/12],["CuS2",-28.68865111,-28.68865111/6],
                            ["ZnS",-33.56936371,-33.56936371/8],["ZnS2",-52.62629248,-52.62629248/12],["Zn8Cu5",-119.03567336,-119.03567336/52],
                            ["Zn35Cu17",-112.33610199,-112.33610199/52],["Cu2GeS3",-123.48186136,-123.48186136/24],["Cu3Ge",-32.92180255,-32.92180255/8],
                            ["GeS",-44.13581527,-44.13581527/8],["GeS2",-272.26921339,-272.26921339/48]]

# CZTS, CZGS & CZSS formation energy for 8 atoms [band gap, VBM, formation energy] From Ground calculation - 8 atoms
CZTS_8 = [1.36, 4.288117, (-38.623184 - 2 * Pure_energies_CZTS[0][2] - Pure_energies_CZTS[1][2] - 4*Pure_energies_CZTS[2][2] - Pure_energies_CZTS[3][2])/8 ,  -38.623184/8]
CZGS_8 = [1.89, 4.021333, (-39.43697516 - 2 * Pure_energies_CZGS[0][2] - Pure_energies_CZGS[1][2] - 4*Pure_energies_CZGS[2][2] - Pure_energies_CZGS[3][2])/8, -39.43697516/8] 
CZTS_64 = [1.18, 3.863666, ((-308.96612677) - 16 * Pure_energies_CZTS[0][2] - 8*Pure_energies_CZTS[1][2] - 32*Pure_energies_CZTS[2][2] - 8*Pure_energies_CZTS[3][2])/64, -308.96612677/64]
CZGS_64 = [1.72, 3.611929, ((-315.53754352) - 16 * Pure_energies_CZGS[0][2] - 8*Pure_energies_CZGS[1][2] - 32*Pure_energies_CZGS[2][2] - 8*Pure_energies_CZGS[3][2])/64, -315.53754352/64]
print("Band gap - VBM - Formation energy - Energy per atom")
print("CZTS \n 8 atoms:" , CZTS_8 , "\n 64 atoms:" , CZTS_64)
print("CZGS \n 8 atoms:" , CZGS_8 , "\n 64 atoms:" , CZGS_64)

#Kesterite and secondary phases formation energies
Formation_CZTS = 8*(-308.96612677/64) - 2 * Pure_energies_CZTS[0][2] - Pure_energies_CZTS[1][2] - 4*Pure_energies_CZTS[2][2] - Pure_energies_CZTS[3][2]
Formation_CZGS = 8*(-315.53754352/64) - 2 * Pure_energies_CZGS[0][2] - Pure_energies_CZGS[1][2] - 4*Pure_energies_CZGS[2][2] - Pure_energies_CZGS[3][2] - 0.27 

# Secondary phases formation energies -> E(Secondary phase) - E(pure elements) for x atoms in the secondary phase (Ex: SnS2 -> 3 atoms)     
Secondaries_formation_CZTS = [["Cu7S4",11*Secondaries_energies_CZTS[0][2] - 7*Pure_energies_CZTS[0][2] - 4*Pure_energies_CZTS[2][2]],
                            ["CuS",2*Secondaries_energies_CZTS[1][2] - Pure_energies_CZTS[0][2] - Pure_energies_CZTS[2][2]],
                            ["CuS2",3*Secondaries_energies_CZTS[2][2] - Pure_energies_CZTS[0][2] - 2*Pure_energies_CZTS[2][2]],
                            ["ZnS",2*Secondaries_energies_CZTS[3][2] - Pure_energies_CZTS[1][2] - Pure_energies_CZTS[2][2]],
                            ["ZnS2",3*Secondaries_energies_CZTS[4][2] - Pure_energies_CZTS[1][2] - 2*Pure_energies_CZTS[2][2]],
                            ["Zn8Cu5",13*Secondaries_energies_CZTS[5][2] - 5*Pure_energies_CZTS[0][2] - 8*Pure_energies_CZTS[1][2]],
                            ["Zn35Cu17",52*Secondaries_energies_CZTS[6][2] - 17*Pure_energies_CZTS[0][2] - 35*Pure_energies_CZTS[1][2]],
                            ["Cu2SnS3",6*Secondaries_energies_CZTS[7][2] - 2*Pure_energies_CZTS[0][2] - 3*Pure_energies_CZTS[2][2] - Pure_energies_CZTS[3][2]],
                            ["CuSn",2*Secondaries_energies_CZTS[8][2] - Pure_energies_CZTS[0][2] - Pure_energies_CZTS[3][2]],
                            ["SnS",2*Secondaries_energies_CZTS[9][2] - Pure_energies_CZTS[2][2] - Pure_energies_CZTS[3][2]],
                            ["SnS2",3*Secondaries_energies_CZTS[10][2] - 2*Pure_energies_CZTS[2][2] - Pure_energies_CZTS[3][2]]]

Secondaries_formation_CZGS = [["Cu7S4",11*Secondaries_energies_CZGS[0][2] - 7*Pure_energies_CZTS[0][2] - 4*Pure_energies_CZTS[2][2]],
                            ["CuS",2*Secondaries_energies_CZGS[1][2] - Pure_energies_CZTS[0][2] - Pure_energies_CZTS[2][2]],
                            ["CuS2",3*Secondaries_energies_CZGS[2][2] - Pure_energies_CZTS[0][2] - 2*Pure_energies_CZTS[2][2]],
                            ["ZnS",2*Secondaries_energies_CZGS[3][2] - Pure_energies_CZTS[1][2] - Pure_energies_CZTS[2][2]],
                            ["ZnS2",3*Secondaries_energies_CZGS[4][2] - Pure_energies_CZTS[1][2] - 2*Pure_energies_CZTS[2][2]],
                            ["Zn8Cu5",13*Secondaries_energies_CZGS[5][2] - 5*Pure_energies_CZTS[0][2] - 8*Pure_energies_CZTS[1][2]],
                            ["Zn35Cu17",52*Secondaries_energies_CZGS[6][2] - 17*Pure_energies_CZTS[0][2] - 35*Pure_energies_CZTS[1][2]],
                            ["Cu2GeS3",6*Secondaries_energies_CZGS[7][2] - 2*Pure_energies_CZGS[0][2] -3*Pure_energies_CZGS[2][2] - Pure_energies_CZGS[3][2]],
                            ["Cu3Ge",4*Secondaries_energies_CZGS[8][2] - 3*Pure_energies_CZGS[0][2] - Pure_energies_CZGS[3][2]],
                            ["GeS",2*Secondaries_energies_CZGS[9][2] - Pure_energies_CZGS[2][2] - Pure_energies_CZGS[3][2]],
                            ["GeS2",3*Secondaries_energies_CZGS[10][2] - 2*Pure_energies_CZGS[2][2] - Pure_energies_CZGS[3][2]]]

# Wexler corrections for Ge secondary phases 
#R. B. Wexler, G. S. Gautam, and E. A. Carter, Exchange-correlation functional challenges in modeling quaternary chalcogenides, Physical Review B 102, 054101 (2020).
Secondaries_formation_CZGS[7][1] -= 0.27
Secondaries_formation_CZGS[8][1] -= 0.27
Secondaries_formation_CZGS[9][1] -= 0.27
Secondaries_formation_CZGS[10][1] -= 0.27

print(Formation_CZTS)
print(Formation_CZGS, '\n')
print("CZTS Formation energies",Secondaries_formation_CZTS)
print(' ')
print("CZGS Formation energies",Secondaries_formation_CZGS)

########################################################################################################################
#                                                                                                                      #
#                                   Phase diagram plot for CZTS                                                        #
#                                                                                                                      #
########################################################################################################################

#_________________________________Plot mu_zn, mucu, mu_sn where mu_s = 0________________________________________________

plt.rcParams['axes.linewidth'] = 1.5
fig_CZTS = plt.figure(figsize=(8, 8))
ax = plt.subplot()

mu_Sn = np.linspace(0, Formation_CZTS, nbr) # Table of 100 pts from 0 to CZTS formation energy 
mu_Zn = np.linspace(0, -5, nbr)
mu_Cu_max = mu_cu_init 
mu_zero = np.linspace(0, 0, nbr)
result_limit = np.zeros(nbr)                                            # Line of the phase diagram limit 
result_secondaries = np.zeros((len(Secondaries_formation_CZTS),nbr))    # Table of secondary phases and calculation points to draw the line
stoichiometry = np.zeros((2,nbr))

def limit(mu_Cu_max): # Function returning the phase diagram limit for mu_S = 0 and for a given value of mu_Cu
    for i in range(nbr):
        result_limit[i] = (Formation_CZTS - mu_Cu_max) - mu_Zn[i]
    ax.plot(mu_Zn,result_limit,'black',label='$\mu_S$ = 0 eV')
    ax.fill_between(mu_Zn,0,result_limit,color="lightgrey", alpha=0.5)
    return ax

def secondaries(mu_Cu_max): # Function calculating the line to be drawn of each secondary phase
    # Calculation of the points in the line
    for i in range(nbr):
        result_secondaries[0][i] = Formation_CZTS - mu_Zn[i] - Secondaries_formation_CZTS[0][1] + 5*mu_Cu_max
        result_secondaries[1][i] = Formation_CZTS - mu_Zn[i] - 4*Secondaries_formation_CZTS[1][1] + 2*mu_Cu_max
        result_secondaries[2][i] = Formation_CZTS - mu_Zn[i] - 2*Secondaries_formation_CZTS[2][1]
        result_secondaries[3][i] = Formation_CZTS + 3*mu_Zn[i] - 4*Secondaries_formation_CZTS[3][1] - 2*mu_Cu_max
        result_secondaries[4][i] = Formation_CZTS + mu_Zn[i] - 2*Secondaries_formation_CZTS[4][1] - 2*mu_Cu_max
        result_secondaries[7][i] = -3*Formation_CZTS + 3*mu_Zn[i] + 4*Secondaries_formation_CZTS[7][1] - 2*mu_Cu_max
        result_secondaries[9][i] = (-Formation_CZTS + mu_Zn[i] + 4*Secondaries_formation_CZTS[9][1] + 2*mu_Cu_max)/3
        result_secondaries[10][i] = -Formation_CZTS + mu_Zn[i] + 2*Secondaries_formation_CZTS[10][1] +2*mu_Cu_max

        # Stochiometric lines whose intersection give the usual kesterite stochiometry
        stoichiometry[0][i] = mu_Zn[i]*1.1
        stoichiometry[1][i] = (mu_Zn[i]*3*mu_Cu_max)/(4*mu_Zn[i] - 3*mu_Cu_max)

    idx = np.argwhere(np.diff(np.sign(stoichiometry[1] - stoichiometry[0]))).flatten()
    mu_Zn_stock = mu_Zn[mu_Zn < 1.5/2*mu_Cu_max] #keep only value above 1.5mu_Cu for stoichiometry
    ndelete = stoichiometry[1][len(mu_Zn)-len(mu_Zn_stock):]

    # Plot of the secondary phases lines
    ax.plot(mu_Zn,result_secondaries[0], 'seagreen'  ,label='$Cu_7S_4$')
    ax.plot(mu_Zn,result_secondaries[1], 'indianred'  , label='$CuS$')
    ax.plot(mu_Zn,result_secondaries[2], 'navy'      , label='$CuS_2$')
    ax.plot(mu_Zn,result_secondaries[3], 'sandybrown', label='$ZnS$')
    ax.plot(mu_Zn,result_secondaries[4], 'steelblue' , label='$ZnS_2$')
    ax.plot(mu_Zn, result_secondaries[7], 'chocolate', label='$Cu_2SnS_3$')
    ax.plot(mu_Zn, result_secondaries[9], 'mediumturquoise', label='$SnS$')
    ax.plot(mu_Zn, result_secondaries[10], 'gold', label='$SnS_2$')

    higherline = np.zeros(nbr)
    lowerline = np.zeros(nbr)
    lowerline = np.amax([result_secondaries[0],result_secondaries[1],result_secondaries[2],result_secondaries[3],result_secondaries[4],result_limit],axis=0)
    higherline = np.amin([result_secondaries[7],result_secondaries[9],result_secondaries[10]],axis=0)
    ax.fill_between(mu_Zn, lowerline, higherline,facecolor='darkgreen', alpha = 0.25, where=(higherline > lowerline))    

    return ax

def phase_diagram(mu_Cu_max): # Call to previous function
    limit(mu_Cu_max)
    secondaries(mu_Cu_max)
    if gif == True:
        axCu = plt.axes([0.15, 0.05, 0.75, 0.03], facecolor='orange')
        cu_slide = Slider(axCu,"$\mu_{Cu}$", Formation_CZTS/4, 0, valinit=mu_Cu_max, valstep=0.001,color='lightgrey')
        cu_slide.label.set_size(15)
        return fig_CZTS

phase_diagram(mu_Cu_max)

#__________________________________________Slider for mu_Cu value_______________________________________________________

axCu = plt.axes([0.15, 0.05, 0.75, 0.03], facecolor='orange')
cu_slide = Slider(axCu,"$\mu_{Cu}$", -1, 0, valinit=mu_Cu_max, valstep=0.001,color='lightgrey')
cu_slide.label.set_size(15)

def update(val): # Update phase diagram for various mu_Cu values
    mu_Cu_max = cu_slide.val
    ax.cla()
    phase_diagram(mu_Cu_max)
    CZTS_plot_look()
    if gif == True:
        return fig_CZTS

cu_slide.on_changed(update)

#________________________________________________Plot look______________________________________________________________
def CZTS_plot_look(): # Phase diagram look
    ax.legend(loc='upper right',prop={'size': 12}) #title="$Cu_2ZnSnS_4$",
    ax.set_xlabel("$\mu_{Zn}$ [eV]",fontdict = fontlabel , labelpad=5)
    ax.set_ylabel("$\mu_{Sn}$ [eV]",fontdict = fontlabel , labelpad=5)
    ax.set_xlim(0,-4.5)
    ax.set_ylim(0,-4.5)
    ax.set_xticks([0,-0.5,-1,-1.5,-2,-2.5,-3,-3.5,-4,-4.5])
    ax.set_yticks([0,-0.5,-1,-1.5,-2,-2.5,-3,-3.5,-4,-4.5])
    ax.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)
    fig_CZTS.subplots_adjust(left=0.15,bottom=0.19,right=0.96,top=0.97)
    if gif == True:
        return fig_CZTS

CZTS_plot_look()

#_______________________________________________GIF generator____________________________________________________________
if gif: 
    valupdate = np.linspace(0,100,50)
    valupdate = -valupdate/100
    camera = Camera(fig_CZTS)
    for i in valupdate:
        phase_diagram(i)
        camera.snap()
    animation = camera.animate() # animation ready
    animation.save(path + "CZTS.gif") 

if Gecompound:
    ########################################################################################################################
    #                                                                                                                      #
    #                                   Phase diagram plot for CZGS                                                        #
    #                                                                                                                      #
    ########################################################################################################################

    #print('Formation energy CZGS',Formation_CZGS)

    #_________________________________Plot mu_zn, mucu, mu_sn where mu_s = 0________________________________________________

    fig_CZGS = plt.figure(figsize=(8, 8))
    ax_CZGS = plt.subplot()

    mu_Sn_CZGS = np.linspace(0, Formation_CZGS, nbr)
    mu_Zn_CZGS = np.linspace(0, -7, nbr)
    mu_Cu_max_CZGS = mu_cu_init
    result_limit_CZGS = np.zeros(nbr)
    result_secondaries_CZGS = np.zeros((len(Secondaries_formation_CZGS),nbr))
    stoichiometry_CZGS = np.zeros((2,nbr))

    def limit_CZGS(mu_Cu_max):
        for i in range(nbr):
            result_limit_CZGS[i] = (Formation_CZGS - mu_Cu_max) - mu_Zn_CZGS[i]

        ax_CZGS.plot(mu_Zn_CZGS,result_limit_CZGS,'black',label='$\mu_S$ = 0 eV')
        ax_CZGS.fill_between(mu_Zn_CZGS,0,result_limit_CZGS,color="lightgrey", alpha=0.5)

    def secondaries_CZGS(mu_Cu_max):
        for i in range(nbr):
            result_secondaries_CZGS[0][i] = Formation_CZGS - mu_Zn_CZGS[i] - Secondaries_formation_CZGS[0][1] + 5*mu_Cu_max
            result_secondaries_CZGS[1][i] = Formation_CZGS - mu_Zn_CZGS[i] - 4*Secondaries_formation_CZGS[1][1] + 2*mu_Cu_max
            result_secondaries_CZGS[2][i] = Formation_CZGS - mu_Zn_CZGS[i] - 2*Secondaries_formation_CZGS[2][1]
            result_secondaries_CZGS[3][i] = Formation_CZGS + 3*mu_Zn_CZGS[i] - 4*Secondaries_formation_CZGS[3][1] - 2*mu_Cu_max
            result_secondaries_CZGS[4][i] = Formation_CZGS + mu_Zn_CZGS[i] - 2*Secondaries_formation_CZGS[4][1] - 2*mu_Cu_max
            result_secondaries_CZGS[7][i] = -3*Formation_CZGS + 3*mu_Zn_CZGS[i] + 4*Secondaries_formation_CZGS[7][1] - 2*mu_Cu_max
            result_secondaries_CZGS[9][i] = (-Formation_CZGS + mu_Zn_CZGS[i] + 4*Secondaries_formation_CZGS[9][1] + 2*mu_Cu_max)/3
            result_secondaries_CZGS[10][i] = -Formation_CZGS + mu_Zn_CZGS[i] + 2*Secondaries_formation_CZGS[10][1] + 2*mu_Cu_max

            stoichiometry_CZGS[0][i] = mu_Zn_CZGS[i]*1.1
            stoichiometry_CZGS[1][i] = (mu_Zn_CZGS[i]*3*mu_Cu_max)/(4*mu_Zn_CZGS[i] - 3*mu_Cu_max)

        #stoch_CZGS = np.linalg.solve([[1, -1 * 1.1], [1, 1]], [0, mu_Cu_max*(2*0.75)])
        idx = np.argwhere(np.diff(np.sign(stoichiometry_CZGS[1] - stoichiometry_CZGS[0]))).flatten()
        mu_Zn_CZGS_stock = mu_Zn_CZGS[mu_Zn_CZGS < 1.5/2*mu_Cu_max] #keep only value above 1.5mu_Cu for stoichiometry -> ccheck eq in notebook
        ndelete = stoichiometry_CZGS[1][len(mu_Zn_CZGS)-len(mu_Zn_CZGS_stock):]

        ax_CZGS.plot(mu_Zn_CZGS,result_secondaries_CZGS[0], 'seagreen'  ,label='$Cu_7S_4$')
        ax_CZGS.plot(mu_Zn_CZGS,result_secondaries_CZGS[1], 'indianred'  , label='$CuS$')
        ax_CZGS.plot(mu_Zn_CZGS,result_secondaries_CZGS[2], 'navy'      , label='$CuS_2$')
        ax_CZGS.plot(mu_Zn_CZGS,result_secondaries_CZGS[3], 'sandybrown', label='$ZnS$')
        ax_CZGS.plot(mu_Zn_CZGS,result_secondaries_CZGS[4], 'steelblue' , label='$ZnS_2$')
        ax_CZGS.plot(mu_Zn_CZGS, result_secondaries_CZGS[7], 'chocolate', label='$Cu_{2}GeS_3$')
        ax_CZGS.plot(mu_Zn_CZGS, result_secondaries_CZGS[9], 'mediumturquoise', label='$GeS$')
        ax_CZGS.plot(mu_Zn_CZGS, result_secondaries_CZGS[10], 'gold', label='$GeS_2$')
        
        idx = idx[2:]
        
        higherline = np.zeros(nbr)
        lowerline = np.zeros(nbr)
        lowerline = np.amax([result_secondaries_CZGS[0],result_secondaries_CZGS[1],result_secondaries_CZGS[2],result_secondaries_CZGS[3],result_secondaries_CZGS[4],result_limit_CZGS],axis=0)
        higherline = np.amin([result_secondaries_CZGS[7],result_secondaries_CZGS[9],result_secondaries_CZGS[10]],axis=0)
        ax_CZGS.fill_between(mu_Zn_CZGS, lowerline, higherline,facecolor='darkgreen', alpha = 0.25, where=(higherline > lowerline))

    def phase_diagram_CZGS(mu_Cu_max_CZGS):
        if gif == True:
            axCu_CZGS = plt.axes([0.15, 0.05, 0.75, 0.03], facecolor='orange')
            cu_slide_CZGS = Slider(axCu_CZGS,"$\mu_{Cu}$", Formation_CZGS, 0, valinit=mu_Cu_max_CZGS, valstep=0.001,color='lightgrey')
            cu_slide_CZGS.label.set_size(15)
        limit_CZGS(mu_Cu_max_CZGS)
        secondaries_CZGS(mu_Cu_max_CZGS)

    phase_diagram_CZGS(mu_Cu_max_CZGS)

    #__________________________________________Slider for mu_Cu value_______________________________________________________
    axCu_CZGS = plt.axes([0.15, 0.05, 0.75, 0.03], facecolor='orange')
    cu_slide_CZGS = Slider(axCu_CZGS,"$\mu_{Cu}$", -1, 0, valinit=mu_Cu_max_CZGS, valstep=0.001,color='lightgrey')
    cu_slide_CZGS.label.set_size(15)

    def update_CZGS(val):
        mu_Cu_max_CZGS = cu_slide_CZGS.val
        ax_CZGS.cla()
        phase_diagram_CZGS(mu_Cu_max_CZGS)
        CZGS_plot_look()

    cu_slide_CZGS.on_changed(update_CZGS)

    #________________________________________________Plot look______________________________________________________________
    def CZGS_plot_look():
        ax_CZGS.legend(loc='upper right',prop={'size': 12}) #title="$Cu_2ZnGeS_4$"
        ax_CZGS.set_xlabel("$\mu_{Zn}$ [eV]",fontdict = fontlabel , labelpad=5)
        ax_CZGS.set_ylabel("$\mu_{Ge}$ [eV]",fontdict = fontlabel , labelpad=5)
        ax_CZGS.set_xlim(0,-4.5)
        ax_CZGS.set_ylim(0,-4.5)
        ax_CZGS.set_xticks([0,-0.5,-1,-1.5,-2,-2.5,-3,-3.5,-4,-4.5])
        ax_CZGS.set_yticks([0,-0.5,-1,-1.5,-2,-2.5,-3,-3.5,-4,-4.5])
        ax_CZGS.tick_params(which='both',length=10,width=1.5,direction='in',pad=5,top=True,right=True,grid_linewidth=0.05,grid_alpha=0.5,labelsize=14)
        fig_CZGS.subplots_adjust(left=0.15,bottom=0.19,right=0.96,top=0.97)

    CZGS_plot_look()
    
    #_______________________________________________GIF generator____________________________________________________________
    if gif == True: 
        camera = Camera(fig_CZGS)
        for i in valupdate:
            phase_diagram_CZGS(i)
            camera.snap()
        animation = camera.animate() # animation ready
        animation.save(path + "CZGS.gif") 

if saveflag:
    fig_CZTS.savefig(path + "CZTS.pdf",facecolor='none',bbox_inches='tight')
    if Gecompound:
        fig_CZGS.savefig(path + "CZGS.pdf",facecolor='none',bbox_inches='tight')

plt.show()
