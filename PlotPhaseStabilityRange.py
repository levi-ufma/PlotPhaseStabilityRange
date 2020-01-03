"""
==========================
Plot phase stability range
==========================
Makes a horizontal bar plot showing the stability range of each phase in an 
open chemical system as a function chemical potential, temperature or voltage.
Created on Sat Jan  05 05:45:00 2019
@authors: Jherfson Castro, Rodolpho Mouta
"""
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from thermodynamics import Thermodynamics as mu2t
from PhaseDiagramOpen import PhaseDiagramOpenAnalyzer
#import PhaseAndPotential
#from mu_to_temp import mu_to_temperature as mu2t

######## Input data ########
# Phases, initial and final chemical potential at which they are stable.
# Bar colors, label colors. Element to which the system is open (Li, Na, O).
#pd = PhaseAndPotential.chempot_range_of_each_phase
element = PhaseDiagramOpenAnalyzer(system=["Li", "Ca", "O"], open_element="O")
pd = element.get_phase_diagram_data()

Phase = list(pd)

mu_i = []
mu_f = []
for potential in pd:
    mu_i.append(pd[potential][0])
    mu_f.append(pd[potential][1])

# print("mu_i: ",mu_i)

# últim  valor do potencial e convertido para temperatura


# LastPotential = mu_f[len(mu_f)-1]

#Phase = ['SnO_2','WO_3','MnO_2','Mn_2O_3','Mn_3O_4','MnWO_4']
#mu_i = [-4.93552791875,-4.93552791875,-4.93552791875,-5.658827941874999,-6.467741869375012,-5.886984983124987]
#mu_f = [-8.017093547499998,-7.554040039999997,-5.658827941874999,-6.467741869375012,-7.376891950000004,-8.030107460000005]
#BarColor = ['orchid','silver','lightblue','orange','indigo','lightgreen', 'navy', 'green', 'teal', 'darkorange','lime']

BarColor = ['lightblue']*len(Phase)

OpenTo = 'O'

if OpenTo == 'O':
    LastPotential = mu2t().print_temperature_corresponding_to_mu_equals(
        mu_f[len(mu_f)-1])['T_Celsius']
else:
    LastPotential = mu_f[len(mu_f)-1]


######## Main plot parameters ########
# Whether or not to convert the chemical potential to temperature (when open
# to oxygen) or to voltage against Li/Li+ (when open to Li) or Na/Na+ (when
# open to Na). The options are 'None', 'T_C', 'T_K', 'V_Li', 'V_Na'.
ConvertTo = 'T_C'
# Range of the plotted quantity. If ConvertTo is set no 'None', the
# quantity is the chemical potential. If it is set to 'T_C' or 'T_K', the
# quantity is the temperature in °C or K, respectively. If it is set to
# 'V_Li' or 'V_Na' the quantity is V vs. Li/Li+ or V vs. Na/Na+, respectively.
Xmin, Xmax = Xlim = 0, LastPotential  # 1350
# Bar height, ranging from 0 to 1 [the spacing between bars
# will be (1 - BarHeight)]; phase label font.
BarHeight = 0.7
PhaseLabelSize = 14
# Y axis range; minor and major ticks spacing; tick label size.
Ymin, Ymax = Ylim = BarHeight-1, len(Phase)-0

# MinorTickSpacing,MajorTickSpacing = 0.5, 2#50, 200


MinorTickSpacing, MajorTickSpacing = 50, 200

TickFontSize = 11
# Axes labels size and their distance from borders.
AxisFontSize = 16

LabelShift = 10
# Grid transparency in %
Transparency = 75


######## Processing of input data ########
# Set Y label
YLabel = 'Stable phases'
# Set X axis label based on the open element and the conversion selected by the user.
if OpenTo == 'O':
    if ConvertTo == 'T_C':
        XLabel = 'Temperature (°C)'
    elif ConvertTo == 'T_K':
        XLabel = 'Temperature (K)'
    elif ConvertTo == 'None':
        XLabel = r'-$\Delta$' + r'$\mu_O$ (eV)'
    else:
        XLabel = r'$\mu_O$ (eV)'
        ConvertTo = 'None'
        print('!!! ERROR !!! \nSelected conversion of oxygen chemical potential not allowed. '
              'The X axis quantity will be kept as chemical potential. '
              "Check 'ConvertTo' variable. ")
elif OpenTo == 'Li':
    if ConvertTo == 'V_Li':
        XLabel = 'V vs. Li/Li$^+$ (V)'
    elif ConvertTo == 'None':
        XLabel = r'$\mu_{Li}$ vs. Li° (eV)'
    else:
        XLabel = r'$\mu_{Li}$ vs. Li° (eV)'
        ConvertTo = 'None'
        print('!!! ERROR !!! \nSelected conversion of lithium chemical potential not allowed. '
              'The X axis quantity will be kept as chemical potential. '
              "Check 'ConvertTo' variable. ")
elif OpenTo == 'Na':
    if ConvertTo == 'V_Na':
        XLabel = 'V vs. Na/Na$^+$ (V)'
    elif ConvertTo == 'None':
        XLabel = r'$\mu_{Na}$ vs. Na° (eV)'
    else:
        XLabel = r'$\mu_{Na}$ vs. Na° (eV)'
        ConvertTo = 'None'
        print('!!! ERROR !!! \nSelected conversion of sodium chemical potential not allowed. '
              'The X axis quantity will be kept as chemical potential. '
              "Check 'ConvertTo' variable.")
else:
    print('!!! ERROR !!!', "\n'{}'".format(OpenTo), 'is not a valid option to '
          "describe to which element the system is open. The only valid "
          "options currently are 'Li', 'Na' and 'O'. Ckeck 'OpenTo' variable.")
    if ConvertTo == 'T_C':
        XLabel = 'Temperature (°C)'
        OpenTo = 'O'
    elif ConvertTo == 'T_K':
        XLabel = 'Temperature (K)'
        OpenTo = 'O'
    elif ConvertTo == 'V_Na':
        XLabel = 'V vs. Na/Na$^+$'
        OpenTo = 'Na'
    elif ConvertTo == 'V_Li':
        XLabel = 'V vs. Li/Li$^+$'
        OpenTo = 'Li'
    elif ConvertTo == 'None':
        XLabel = r'$\mu$' + r'$\mu$ (eV)'
    else:
        XLabel = r'$\mu$ (eV)'
        ConvertTo = 'None'
    if ConvertTo == 'None':
        print("Thus, the x axis variable will be assumed to be an arbitrary "
              "chemical potential.")
    else:
        print("However, 'ConvertTo' variable is set to", "'{}'".format(ConvertTo), 'which '
              'is associated to', '{}.'.format(
                  OpenTo), " Thus, it was assumed that the "
              'system is in fact open to', '{},'.format(
                  OpenTo), "so that the 'OpenTo' "
              "variable was overridden and set to this element.")
# Convert the the data from chemical potential to temperature or voltage,
# if necessary (i.e., if 'ConvertTo' is not set to 'None').
if ConvertTo == 'T_C':
    Xi = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Celsius']
          for n in mu_i]
    Xf = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Celsius']
          for n in mu_f]
elif ConvertTo == 'T_K':
    Xi = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Kelvin']
          for n in mu_i]
    Xf = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Kelvin']
          for n in mu_f]
elif ConvertTo == 'V_Li':
    # -1.908 is the chemical potential of Li at 0 K.
    Xi = [-1.908 - n for n in mu_i]
    # -1.908 is the chemical potential of Li at 0 K.
    Xf = [-1.908 - n for n in mu_f]
elif ConvertTo == 'V_Na':
    # -1.313 is the chemical potential of Na at 0 K.
    Xi = [-1.313 - n for n in mu_i]
    # -1.313 is the chemical potential of Na at 0 K.
    Xf = [-1.313 - n for n in mu_f]
elif ConvertTo == 'None':
    if OpenTo == 'Li':
        Xi = [n - (-1.908) for n in mu_i]
        Xf = [n - (-1.908) for n in mu_f]
    elif OpenTo == 'Na':
        Xi = [n - (-1.313) for n in mu_i]
        Xf = [n - (-1.313) for n in mu_f]
    else:
        Xi = [mu_i[0] - n for n in mu_i]
        print(Xi)
        Xf = [mu_i[0] - n for n in mu_f]
        print(Xf)

# Initialize label color. Trim initial and final X values, based on the
# X range. Also, if the bar color is too dark, change the label
# color to white instead of black. Put the phase labels in boldface.
LabelColor = ['black']*len(Phase)
for n in range(0, len(Phase)):
    if Xi[n] < Xmin:
        Xi[n] = Xmin
    if Xf[n] > Xmax:
        Xf[n] = Xmax
    # if BarColor[n] in ['black','Black','brown','Brown','navy','Navy','blue','Blue','green','Green', 'red', 'Red', 'green', 'Green', 'lime', 'Lime', 'indigo', 'Indigo', 'navy', 'Navy', 'darkorange', 'Darkorange']:
        #LabelColor[n] = 'white'
    Phase[n] = '${'+Phase[n]+'}$'
#Xf[1] = 563.45
# Calculate range and center of X axis quantity for each phase. Generate y positions.
Xr = (np.array(Xf) - np.array(Xi)).tolist()
Xc = ((np.array(Xf) + np.array(Xi))/2).tolist()
BarYPos = np.arange(len(Phase))
# BarYPos[3] = BarYPos[2]
# BarYPos[4] = BarYPos[2]
# BarYPos[5] = BarYPos[1]


######## Plot generation ########
fig, ax = plt.subplots(figsize=(6, 2), dpi=115)
# Set plot limits.
ax.set(xlim=Xlim, ylim=Ylim, autoscale_on=False)
# Generate bars.
ax.barh(y=BarYPos, left=Xi, width=Xr,
        height=BarHeight, align='edge', color=BarColor)
# Generate phase labels.
LabelProps = {'horizontalalignment': 'center', 'verticalalignment': 'center',
              'fontsize': PhaseLabelSize, 'fontweight': 'bold'}
for n in range(0, len(Phase)):
    ax.text(x=Xc[n], y=BarYPos[n]+BarHeight/2,
            s=Phase[n], color=LabelColor[n], **LabelProps)
# Generate axes labels.
ax.set_xlabel(XLabel, fontsize=AxisFontSize, labelpad=LabelShift)
ax.set_ylabel(YLabel, fontsize=AxisFontSize, labelpad=LabelShift)
# Set major and minor ticks spacing, direction, and respective label size
ax.xaxis.set_major_locator(plt.MultipleLocator(MajorTickSpacing))
ax.xaxis.set_minor_locator(plt.MultipleLocator(MinorTickSpacing))
ax.tick_params(axis='both', which='both', direction='in',
               labelsize=TickFontSize, width=0.8)
ax.tick_params(axis='both', which='minor', length=3)
ax.tick_params(axis='both', which='major', length=5)
ax.set_yticks([])
# Set vertical grid lines. Choose if they follow only major ticks,
# minor ticks or both. Set grid transparency.
ax.grid(axis='x', which='minor', alpha=1-Transparency/100)

#plt.savefig("LiCaO-Li+", dpi=1200)
plt.show()
