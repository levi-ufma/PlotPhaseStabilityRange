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
import numpy as np
from mu_to_temp2_2 import mu_to_temperature as mu2t

######## Input data ########
# Phases, initial and final chemical potential at which they are stable. 
# Bar colors, label colors. Element to which the system is open (Li, Na, O).
Phase = ['SnO_2','WO_3','MnO_2','Mn_2O_3','Mn_3O_4','MnWO_4']
mu_i = [-4.936,-4.936,-4.936,-5.659,-6.468,-5.887]
mu_f = [-8.017,-7.554,-5.659,-6.468,-7.377,-8.030]
BarColor = ['orchid','silver','lightblue','orange','tomato','lightgreen']
OpenTo = 'O'


######## Main plot parameters ########
# Whether or not to convert the chemical potential to temperature (when open 
# to oxygen) or to voltage against Li/Li+ (when open to Li) or Na/Na+ (when
# open to Na). The options are 'None', 'T_C', 'T_K', 'V_Li', 'V_Na'.
ConvertTo = 'T_C'
# Range of the plotted quantity. If ConvertTo is set no 'None', the
# quantity is the chemical potential. If it is set to 'T_C' or 'T_K', the 
# quantity is the temperature in °C or K, respectively. If it is set to 
# 'V_Li' or 'V_Na' the quantity is V vs. Li/Li+ or V vs. Na/Na+, respectively. 
Xmin,Xmax = Xlim = 0, 1350
# Bar height, ranging from 0 to 1 [the spacing between bars 
# will be (1 - BarHeight)]; phase label font.
BarHeight = 0.7
PhaseLabelSize = 14
# Y axis range; minor and major ticks spacing; tick label size.
Ymin,Ymax = Ylim = BarHeight-1, len(Phase)
MinorTickSpacing,MajorTickSpacing = 50, 150
TickFontSize = 11
#Axes labels size and their distance from borders.
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
        XLabel = '$\mu_O$ (eV)'
    else:
        XLabel = '$\mu_O$ (eV)'
        ConvertTo = 'None'
        print('!!! ERROR !!! \nSelected conversion of oxygen chemical potential not allowed. '
              'The X axis quantity will be kept as chemical potential. '
              "Check 'ConvertTo' variable. ")
elif OpenTo == 'Li':
    if ConvertTo == 'V_Li':
        XLabel = 'V vs. Li/Li$^+$ (V)'
    elif ConvertTo == 'None':
        XLabel = '$\mu_{Li}$ vs. Li° (eV)'
    else:
        XLabel = '$\mu_{Li}$ vs. Li° (eV)'
        ConvertTo = 'None'
        print('!!! ERROR !!! \nSelected conversion of lithium chemical potential not allowed. '
              'The X axis quantity will be kept as chemical potential. '
              "Check 'ConvertTo' variable. ")
elif OpenTo == 'Na':
    if ConvertTo == 'V_Na':
        XLabel = 'V vs. Na/Na$^+$ (V)'
    elif ConvertTo == 'None':
        XLabel = '$\mu_{Na}$ vs. Na° (eV)'
    else:
        XLabel = '$\mu_{Na}$ vs. Na° (eV)'
        ConvertTo = 'None'
        print('!!! ERROR !!! \nSelected conversion of sodium chemical potential not allowed. '
              'The X axis quantity will be kept as chemical potential. '
              "Check 'ConvertTo' variable.")
else:
    print('!!! ERROR !!!',"\n'{}'".format(OpenTo),'is not a valid option to '
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
        XLabel = '$\mu$ (eV)'
    else:
        XLabel = '$\mu$ (eV)'
        ConvertTo = 'None'
    if ConvertTo == 'None':
        print("Thus, the x axis variable will be assumed to be an arbitrary "
              "chemical potential.")
    else:
        print("However, 'ConvertTo' variable is set to","'{}'".format(ConvertTo),'which '
              'is associated to','{}.'.format(OpenTo)," Thus, it was assumed that the "
              'system is in fact open to','{},'.format(OpenTo),"so that the 'OpenTo' "
              "variable was overridden and set to this element.")
# Convert the the data from chemical potential to temperature or voltage, 
# if necessary (i.e., if 'ConvertTo' is not set to 'None').
if ConvertTo == 'T_C':
    Xi = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Celsius'] for n in mu_i]
    Xf = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Celsius'] for n in mu_f]
elif ConvertTo == 'T_K':
    Xi = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Kelvin'] for n in mu_i]    
    Xf = [mu2t().print_temperature_corresponding_to_mu_equals(n)['T_Kelvin'] for n in mu_f]    
elif ConvertTo == 'V_Li':
    Xi = [-1.908 - n for n in mu_i] # -1.908 is the chemical potential of Li at 0 K.
    Xf = [-1.908 - n for n in mu_f] # -1.908 is the chemical potential of Li at 0 K.
elif ConvertTo == 'V_Na':
    Xi = [-1.313 - n for n in mu_i] # -1.313 is the chemical potential of Na at 0 K.
    Xf = [-1.313 - n for n in mu_f] # -1.313 is the chemical potential of Na at 0 K.
elif ConvertTo == 'None':
    if OpenTo == 'Li':
        Xi = [n - (-1.908) for n in mu_i]
        Xf = [n - (-1.908) for n in mu_f]
    elif OpenTo == 'Na':
        Xi = [n - (-1.313) for n in mu_i]
        Xf = [n - (-1.313) for n in mu_f]
    else:
        Xi = mu_i
        Xf = mu_f
# Initialize label color. Trim initial and final X values, based on the
# X range. Also, if the bar color is too dark, change the label
# color to white instead of black. Put the phase labels in boldface.
LabelColor = ['black']*len(Phase)
for n in range(0,len(Phase)):
    if Xi[n] < Xmin:
        Xi[n] = Xmin
    if Xf[n] > Xmax:
        Xf[n] = Xmax
    if BarColor[n] in ['black','Black','brown','Brown','navy','Navy','blue','Blue','green','Green']:
        LabelColor[n] = 'white'
    Phase[n]='$\mathbf{'+Phase[n]+'}$'
# Calculate range and center of X axis quantity for each phase. Generate y positions.
Xr=(np.array(Xf) - np.array(Xi)).tolist()
Xc=((np.array(Xf) + np.array(Xi))/2).tolist()
BarYPos=np.arange(len(Phase))


######## Plot generation ########
fig, ax = plt.subplots()
# Set plot limits.
ax.set(xlim=Xlim, ylim=Ylim, autoscale_on=False)
# Generate bars.
ax.barh(y=BarYPos, left=Xi, width=Xr, height=BarHeight, align='edge', color=BarColor)
# Generate phase labels.
LabelProps = {'horizontalalignment': 'center', 'verticalalignment': 'center', 'fontsize': PhaseLabelSize, 'fontweight':'bold'}
for n in range(0,len(Phase)):
    ax.text(x=Xc[n], y=BarYPos[n]+BarHeight/2, s=Phase[n], color=LabelColor[n], **LabelProps)
# Generate axes labels.
ax.set_xlabel(XLabel, fontsize=AxisFontSize, labelpad = LabelShift)
ax.set_ylabel(YLabel, fontsize=AxisFontSize, labelpad = LabelShift)
# Set major and minor ticks spacing, direction, and respective label size
ax.xaxis.set_major_locator(plt.MultipleLocator(MajorTickSpacing))
ax.xaxis.set_minor_locator(plt.MultipleLocator(MinorTickSpacing))
ax.tick_params(axis='both',which='both', direction='in', labelsize=TickFontSize)
ax.set_yticks([])
# Set vertical grid lines. Choose if they follow only major ticks,
# minor ticks or both. Set grid transparency.
ax.grid(axis='x', which='minor', alpha=1-Transparency/100)

plt.show()
