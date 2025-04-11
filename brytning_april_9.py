import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as stat

"""
Brytningsindex med hjälp av laseravståndsmätare
"""

#vätskedjupet i meter
x = np.array([2.85,5.34,7.75,10.22,12.62,14.9])*10**(-2)
dx = 0.001

#avstånd mätt med lasern
y = np.array([0.945,0.954,0.963,0.970,0.982,0.986])
dy = 0.003

#plt.errorbar(x, y, xerr=dx, yerr=dy, fmt='r ')
#plt.plot(x,y,'r.')

slope, intercept, r, p, se = stat.linregress(x, y)
q = np.linspace(0.8*min(x),1.2*max(x),100)
#plt.plot(q,slope*q+intercept, label = f'y = {round(slope,3)}x + {round(intercept,3)}')

#plt.xlim(0.9*min(x),1.05*max(x))
#plt.xlabel('x/m')
#plt.ylabel('y/m')
#plt.ylim(0.98*min(y),1.02*max(y))
#plt.legend(loc=4)
#plt.show()

#################
# Beräkningar
#################

H = 0.865+0.0428
dH = 0.003

Δy = 0.935
dΔy = 0.003

θ1 = np.arccos(H/Δy)*180/np.pi #infallsvinkel i grader
dθ1 = (dH/H+dΔy/Δy)/np.tan(θ1)*180/np.pi #osäkerheten i infallsvinkeln

# För beräkningar i Python används radianer
θ1 = np.arccos(H/Δy) #infallsvinkel i radianer

#praktisk parameter för senare beräkningar
p = slope + 1/np.cos(θ1)

#regression
reg = np.polyfit(x,y,deg=1,cov=True)
alpha = reg[0][0]
konst = reg[0][1]
s_alpha = np.sqrt(reg[1][0,0])
s_konst = np.sqrt(reg[1][1,1])

dp = s_alpha+np.pi*np.sin(θ1)/(180*(np.cos(θ1))**2)*dθ1

tal = p**2*np.sqrt(1-4*(np.sin(θ1))**2/p**2)+p**2-2*(np.sin(θ1))**2
nam = np.sqrt(1-4*(np.sin(θ1))**2/p**2)*np.sqrt(1+np.sqrt(1-4*(np.sin(θ1))**2/p**2))
dn2dp = 1/(np.sqrt(2)*p**2)*tal/nam

tal = -np.sqrt(2*p)*np.pi*np.sin(θ1)*np.cos(θ1)
nam = 180*np.sqrt(p+np.sqrt(p**2-4*(np.sin(θ1))**2))*np.sqrt(p**2-4*(np.sin(θ1))**2)
dn2dt = tal/nam

#brytningsindex för vatten
D1 = 1-(2*np.sin(θ1)/p)**2
D2m = 1-np.sqrt(D1)
D2p = 1+np.sqrt(D1)
nw = 1/np.sqrt(2)*p*np.sqrt(D2p) #lösning nw<1 är ofysikalisk

dnw = np.sqrt(dn2dp**2*dp**2+dn2dt**2*dθ1**2)

print(f'Brytningsindexet för vattten är {np.around(nw,3)}  {np.around(dnw,3)}')

x_fit = np.arange(0.02,0.2,0.01)
y_fit = konst+alpha*x_fit

plt.errorbar(x,y,xerr=dx, yerr=dy,fmt='.', label='data', color='black')
plt.xlabel('x [m]')
plt.xlim(0.9*min(x),1.05*max(x))
plt.ylabel('y [m]')
plt.ylim(0.99*min(y),1.01*max(y))

plt.plot(x_fit,y_fit,'orange', label= f'y = {np.around(alpha,3)}x + {np.around(konst,3)}')
plt.legend()
plt.grid()
plt.title('Uppmätt längd med bygglasern mot vattendjupet')
#plt.fill_between(x_fit,y_fit+s_alpha, y_fit-s_alpha, alpha=0.5,color='grey')
#SER UT ATT VARA UNDERLIGT STORA OSÄKERHETER I LINJEN?

plt.show()