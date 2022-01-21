#!/usr/bin/python
import matplotlib.pyplot as plt

import pylab as pl
import sys
import numpy as np

from scipy.linalg import inv
from matplotlib.widgets import Slider, Button, RadioButtons

from matplotlib.pyplot import text


# IMPORT THE CC, C1, C2 directly from CST

# FREQUENCY VECTOR
pi=np.pi

# Best <#Parameters#> so far

EL_10=200.5
EL_20=150
LL0=659.5
L_SAMPLE0=76.5 # Mesure @ CST 6 turns
C_10=6.32
C_20=91.85
C0=3.85
C_TX0=23.58
L_MX0=77.85

# simply in order to keep the aspect ratio... to caculate BS (to be destroyed by fire and sulfure)
Freq=np.ones_like(5000) 

# Define a Cursor
fig, ax1 = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
freq = np.linspace(5, 245, 5000)




def print_s (EL_1, EL_2, LL, L_SAMPLE, C_1, C_2, C, C_TX, L_MX, Freq):
    
    omega=2.0*pi*Freq                                                            # angular momentum variation range
    frq=Freq*1e-6

    # RESONNANCE FREQUENCY
    f0=400.0                                                                     # resonance frequency (in MHz)
    #f0=input('resonance frequency (MHz):')
    f0=f0*1e6                                                                  # resonance frequency in Hz and
    om0=2.0*pi*f0                                                              # equivalent angular momentum in rad/s
    
    # REDOR BOX PARAMETERS
    
    #The conventions are in VIPEC : \COUPLERS-UCCS\vipec-400\H8906_V5.ckt
    
    #L_in=5.915e-07
    L_in=LL*1e-9              # Passe par parametres
    
    
    #C_1=39.33
    C1=C_1*1e-12              # Passe par parametres
    Zc1 = 1.0/(1j*C1*omega)
    
    
    #C_2=69.85
    C2=C_2*1e-12                  # Passe par parametres
    Zc2 = 1.0/(1j*C2*omega)
    
    # Splitter parameters
    #C=4.40
    CHV=C*1e-12                 # Passe par parametres

    Q_redorC=700.0
              # Ne Passe pas par parametres ajustables (variable locale)
   
   # Construction de la Matrice REDOR
   
    Zc0=1.0/(1j*CHV*omega+CHV*omega/Q_redorC)+1j*L_in*omega
    # <#matrix_REDOR#>
    Tredor11=np.ones_like(omega)+Zc0/Zc1
    Tredor12=-Zc0
    Tredor21=-(1.0/Zc1+1.0/Zc2+Zc0/(Zc1*Zc2))
    Tredor22=np.ones_like(omega)+Zc0/Zc2
    

    # REDOR BOX to X channel Resonant transmission line 50

# PARAMETERS OF UNMATCHED (resonannt) TRANSMISSION LINE (T1) (convention "left to right currents")

    Z1=50.0                                                                    # impedance of tlin (in Ohm)
    #Z2=input('impedance of tlin (Ohm):')
    alpha1=1e-4                                                                # attenuation coefficient for tlin, is a teflon based line
    #alpha2=input('alpha2:')
    # Definition of an electrical lengh modulo 2pi (to justify low cable resonance ~ 70MHz)
    # Harmonic cable resonnace at 241 MHz coresponds to a pi but is not linear
    
    DEG=EL_1
    E_rad1=DEG*pi/180.0                                                          # wavelength in tlin (in radians)
    #E_rad2=input('wavelength in tlin (rad):')
    
    beta1=E_rad1*omega/om0                                                   # phase change coefficient for tlin
    gamma1=alpha1+1j*beta1                                                     # propagation coefficient for tlin
    
    
    # <#coax_t1#>
    T1_11=np.cosh(gamma1)
    T1_12=-Z1*np.sinh(gamma1)
    T1_21=-(1.0/Z1)*np.sinh(gamma1)
    T1_22=np.cosh(gamma1)


    # PARAMETERS OF STOP CIrcuit and capt to ground

    # transfer matrix Ttank (convention "left to right currents")
    # PARAMETERS OF STOP INDUCTANCE
    Cs=5.6*1e-12
    Q_cs=50

    YCs=1j*Cs*omega+Cs*omega/Q_cs                          # admittance of stop capacitance

    #Ls=1.0/(Cs*om0**2)# inductance value (in H)
    Ls=29.5*1e-9
    Qs=500                                                  # quality factor of inductance
    #QL=input('quality factor of inductance:')

    YLs=1.0/(1j*Ls*omega)+1.0/(Ls*Qs*omega)                 # admittance of stop inductance

                                                   # Tunning capacitor CHV in pF
    CTX=C_TX*1e-12
    Zccc = 1.0/(1j*CTX*omega)+1.0/(YLs+YCs)

    Cg=0.5                                                  # Capacitor to ground (maybee repalced)
    Cg=Cg*1e-12

    ZcMH = 1.0/(1j*Cg*omega)

                                                             # matching inductance passed as param
    LMX = L_MX*1e-9                                           # matching inductance value (in H)

    ZMX=1j*LMX*omega                                         #  impedance of matching  inductance

    # <#matrix_Xchanel#> including tank circuit

    Ttloc11=np.ones_like(omega)+Zccc/ZcMH
    Ttloc12=-Zccc
    Ttloc21=-(1.0/ZcMH+1.0/ZMX+Zccc/(ZcMH*ZMX))
    Ttloc22=np.ones_like(omega)+Zccc/ZMX


    # PARAMETERS OF INDUCTANCE

    L=L_SAMPLE                                                      # inductance value (in nH)
    #L=input('inductance value (nH):')
    L=L*1e-9                                                   # inductance value in H
    QL=5000                                                 # quality factor of inductance

    YL=1.0/(1j*L*omega)+1.0/(L*QL*omega)                         # admittance of inductance
    #RL=L*omega*QL                                             # equivalent resistance in parallel

    #transfer matrix TL (convention "left to right currents")
    TL11=np.ones_like(omega)
    TL12=-np.ones_like(omega)/YL
    TL21=0*np.ones_like(omega)
    TL22=np.ones_like(omega)


    # PARAMETERS OF TRANSMISSION LINE (TLIN)

    Z2=40.0                                                                      # impedance of tlin (in Ohm)
    #Z2=input('impedance of tlin (Ohm):')
    alpha2=1e-5                                                                # attenuation coefficient for tlin
    #alpha2=input('alpha2:')
    E_rad2=EL_2*pi/180.0                                                          # wavelength in tlin (in radians)
    #E_rad2=input('wavelength in tlin (rad):')

    beta2=E_rad2*omega/om0                                                   # phase change coefficient for tlin
    gamma2=alpha2+1j*beta2                                                     # propagation coefficient for tlin

    # transfer matrix Ttlin (convention "left to right currents")
    Ttlin11=np.cosh(gamma2)
    Ttlin12=-Z2*np.sinh(gamma2)
    Ttlin21=-(1.0/Z2)*np.sinh(gamma2)
    Ttlin22=np.cosh(gamma2)


    # PARAMETERS OF THE PI-CIRCUIT

    CC=2.2                                                        # coupling capacity value (in pF)
    #CC=input('coupling capacity value (pF):')
    CC=CC*1e-12                                                   # coupling capacity value in F
    QCC=300.0                                                     # quality factor of coupling capacity
    #QCC=input('quality factor of coupling capacity:')

    ZCC=1.0/(1j*CC*omega+CC*omega/QCC)                             # impedance of coupling capacity
    RCC=QCC/(CC*omega)                                             # equivalent resistance in parallel

    CTH=7.725                                                      # tuning capacity value (in pF)
    #CTH=input('tuning capacity value (pF):')
    CTH=CTH*1e-12                                                  # tuning capacity value in F
    QCTH=500.0                                                     # quality factor of tuning capacity
    #QCTH=input('quality factor of tuning capacity:')

    ZCTH=1.0/(1j*CTH*omega+CTH*omega/QCTH)                         # impedance of tuning capacity
    RCTH=QCTH/(CTH*omega)                                          # equivalent resistance in parallel

    CMH=46.535                                                     # matching capacity value (in pF)
    #CMH=input('matching capacity value (pF):')
    CMH=CMH*1e-12                                                  # matching capacity value in F
    QCMH=200.0                                                     # quality factor of matching capacity
    #QCMH=input('quality factor of matching capacity:')

    ZCMH=1.0/(1j*CMH*omega+CMH*omega/QCMH)                         # impedance of matching capacity
    RCMH=QCMH/CMH/omega                                            # equivalent resistance in parallel

    #transfer matrix TC (convention "entering currents")
    TC11=1+ZCC/ZCTH
    TC12=-ZCC
    TC21=(1.0/ZCMH+1.0/ZCTH+ZCC/(ZCTH*ZCMH))
    TC22=-(1+ZCC/ZCMH)


    # INTERMEDIATE TRANSFER MATRIX T (T1 * Tredor )
    T011=T1_11*Tredor11+T1_12*Tredor21
    T012=T1_11*Tredor12+T1_12*Tredor22
    T021=T1_21*Tredor11+T1_22*Tredor21
    T022=T1_21*Tredor12+T1_22*Tredor22
    
    # INTERMEDIATE TRANSFER MATRIX T (Ttloc * T0 )
    T11=Ttloc11*T011+Ttloc12*T021
    T12=Ttloc11*T012+Ttloc12*T022
    T21=Ttloc21*T011+Ttloc22*T021
    T22=Ttloc21*T012+Ttloc22*T022
    

    # INTERMEDIATE TRANSFER MATRIX TI (TL * T)
    TI11=TL11*T11+TL12*T21
    TI12=TL11*T12+TL12*T22
    TI21=TL21*T11+TL22*T21
    TI22=TL21*T12+TL22*T22


    # INTERMEDIATE TRANSFER MATRIX TII (Ttlin * TI)
    TII11=Ttlin11*TI11+Ttlin12*TI21
    TII12=Ttlin11*TI12+Ttlin12*TI22
    TII21=Ttlin21*TI11+Ttlin22*TI21
    TII22=Ttlin21*TI12+Ttlin22*TI22


    # FINAL TRANSFER MATRIX TF (TC * TII)
    TF11=TC11*TII11+TC12*TII21
    TF12=TC11*TII12+TC12*TII22
    TF21=TC21*TII11+TC22*TII21
    TF22=TC21*TII12+TC22*TII22
    detTF=TF11*TF22-TF12*TF21


    # IMPEDANCE MATRIX OF THE PROBE
    Z11=-TF22/TF21
    Z12=1.0/TF21
    Z21=-detTF/TF21
    Z22=TF11/TF21



    # SCATTERING MATRIX S OF THE PROBE

    Z01=50.0                                                                 # reference impedances
    #Z01=377                                                                   # at the two ports
    Z02=50.0                                                                     # of the circuit

    DEN=(1.0+Z11/Z01)*(1.0+Z22/Z02)-(Z12*Z21/(Z01*Z02))

    S11=(1.0-Z11/Z01)*(1.0+Z22/Z02)+(Z12*Z21/(Z01*Z02))
    S11=S11/DEN

    S22=(1.0+Z11/Z01)*(1.0-Z22/Z02)+(Z12*Z21/(Z01*Z02))
    S22=-S22/DEN

    S12=Z12/Z01
    S12=-S12/DEN

    S21=Z21/Z02
    S21=S21/DEN
    
    #print S22;
    return S11;


font = {'family' : 'serif',
       'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
}
# Create a new subplot from a grid of 1x1


# FREQUENCY VECTOR
dim=5000                                                                # number of points, this will be changed
#dim=input('number of points :')
start=5                                                                # start frequency (in MHz)
#start=input('start frequency (MHz):')
startHz=start*1e6                                                         # start frequency in Hz
stop=250                                                                # stop frequency (in MHz)
#stop=input('stop frequency (MHz):')
stopHz=stop*1e6                                                           # stop frequency in Hz

freq=np.linspace(startHz,stopHz,dim)                                      # frequency variation range                                                        # angular momentum variation range
frq=freq*1e-6


# <#best#>
#S11  = print_s (LL=610, C_1=60.725, C_2=3.2, C=4.25, C_TX=50, L_MX=9.5, Freq=freq)
S11  = print_s (EL_1=EL_10, EL_2=EL_20, LL=LL0, L_SAMPLE=L_SAMPLE0, C_1=C_10, C_2=C_20, C=C0, C_TX=C_TX0, L_MX=L_MX0, Freq=freq)
S11_ref=20*np.log10(abs(S11))

## We assume that the BS removes 20% of the field in the coil
S11  = print_s (EL_1=EL_10, EL_2=EL_20, LL=LL0, L_SAMPLE=L_SAMPLE0*0.8, C_1=C_10, C_2=C_20, C=C0, C_TX=C_TX0, L_MX=L_MX0, Freq=freq)
S11_bs = 20*np.log10(abs(S11))
#
#start_search=2030
#range_search=200
#
#boundary_search_range=np.argmin(S11_ref[start_search:start_search+range_search])
##print "minimum freqency"
##print frq[boundary_search_range+start_search]
##print "absolute minimum value"
##print S11_ref[boundary_search_range+start_search]
##ref_freq_low=np.argmin(S11_bs[boundary_search_range+start_search+50:boundary_search_range+start_search+range_search])
##ref_freq_low=np.argmin(S11_ref[boundary_search_range+start_search-range_search:boundary_search_range+start_search-50])
#ref_freq_low=frq[boundary_search_range+start_search]
#print "reference freqency low"
#print ref_freq_low
#
#ref_freq_high=np.argmin(S11_ref[boundary_search_range+start_search-range_search:boundary_search_range+start_search-50])
##ref_freq_high=np.argmin(S11_bs[boundary_search_range+start_search+50:boundary_search_range+start_search+range_search])
#print "reference freqency high"
#print frq[ref_freq_high+start_search]
#ref_freq_BS=frq[ref_freq_high+start_search]
#print "Ball-shift for the Sample resonance"
#BALL_SHIFT=ref_freq_BS-ref_freq_low
#print BALL_SHIFT


ax1=plt.subplot(111)
#plt.axis([50, 150, -10, 10])
l, = ax1.plot(frq, S11_ref , color="red", linewidth=2, linestyle="-", label="reference S11")
ax1.plot(frq, S11_bs, color="blue", linewidth=2, linestyle="-", label="ball-shift")
ax1.set_xlabel('Freq in MHz', fontdict=font)
ax1.set_ylabel('Reflection from port X in dB', color='k', fontdict=font)
#ax1.legend(bbox_to_anchor=(0., 1.08, 1., .108), loc=3,ncol=2, mode="expand", borderaxespad=0.)
#ax1.text(frq[boundary_search_range+start_search]+40,S11_ref[boundary_search_range+start_search]-1, 'Ball-shift for the Sample resonance' , horizontalalignment='center', fontsize=14)
#ax1.text(frq[boundary_search_range+start_search]+40,S11_ref[boundary_search_range+start_search]-2.2,'{: f} MHz'.format(BALL_SHIFT), horizontalalignment='center', fontsize=14)
axcolor = 'lightgoldenrodyellow'
axC1 = plt.axes([0.25, 0.1, 0.65, 0.03])
axC2 = plt.axes([0.25, 0.15, 0.65, 0.03])
axC = plt.axes([0.25, 0.2, 0.65, 0.03])
axLL = plt.axes([0.25, 0.25, 0.65, 0.03])
axCTX = plt.axes([0.25, 0.3, 0.65, 0.03])
axLMX = plt.axes([0.25, 0.35, 0.65, 0.03])
axEL1 = plt.axes([0.25, 0.4, 0.65, 0.03])
axEL2 = plt.axes([0.25, 0.45, 0.65, 0.03])

sC1 = Slider(axC1, 'C1', 2, 122.0, valinit=C_10)
sC2 = Slider(axC2, 'C2', 2, 122.0, valinit=C_20)
sC = Slider(axC, 'C', 0.1, 10.0, valinit=C0)
sLL = Slider(axLL, 'LL', 550.0, 680.0, valinit=LL0)
sCTX = Slider(axCTX, 'CTX', 0.1, 45.0, valinit=C_TX0)
sLMX = Slider(axLMX, 'LMX', 1, 80, valinit=L_MX0)
sEL1 = Slider(axEL1, 'EL1', 20, 390, valinit=EL_10)
sEL2 = Slider(axEL2, 'EL2', 100, 180, valinit=EL_20)

def update(val):
    C1 = sC1.val
    C2 = sC2.val
    C = sC.val
    LL = sLL.val
    CTX = sCTX.val
    LMX = sLMX.val
    EL1 = sEL1.val
    EL2 = sEL2.val
    
    l.set_ydata(20*np.log10( print_s (EL_1=EL1, EL_2=EL2, LL=LL, L_SAMPLE=73.5, C_1=C1, C_2=C2, C=C, C_TX=CTX, L_MX=LMX, Freq=freq)) )
    fig.canvas.draw_idle()
sC1.on_changed(update)
sC2.on_changed(update)
sC.on_changed(update)
sLL.on_changed(update)
sCTX.on_changed(update)
sLMX.on_changed(update)
sEL1.on_changed(update)
sEL2.on_changed(update)

plt.show()
