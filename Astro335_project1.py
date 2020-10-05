import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

Time= np.genfromtxt('C:/Users/Sumner/Downloads/Planet_Finder_Data.csv',delimiter=',', skip_header=1, usecols=0)
relative_flux= np.genfromtxt('C:/Users/Sumner/Downloads/Planet_Finder_Data.csv',delimiter=',', skip_header=1, usecols=1)
error = np.genfromtxt('C:/Users/Sumner/Downloads/Planet_Finder_Data.csv',delimiter=',', skip_header=1, usecols=2)

Sx = 0
Sy = 0
Sr = 1

Py = 0

Pr = 0.125
TimeScale = 17
u = 0.6

def ilim( SunRad, dist, u ):
    return 1 - u * (1 - np.sqrt( (SunRad**2 - dist**2) / SunRad**2))

def ratio_sun_vis(x, SunRad):
    return (np.pi * SunRad**2 - x)/(np.pi * SunRad ** 2)

def overlap_A( dist, SunRad, PlanetRad ):
    radSsqr = SunRad ** 2
    radPsqr = PlanetRad ** 2
   
    angle1 = (radSsqr + dist ** 2 - radPsqr) / (2 * SunRad * dist)
    angle2 = (radPsqr + dist ** 2 - radSsqr) / (2 * PlanetRad * dist)

    theta1 = 2.0 * np.arccos(angle1)
    theta2 = 2.0 * np.arccos(angle2)
   
    area1 = (theta2 - np.sin(theta2)) * (radPsqr / 2.0)
    area2 = (theta1 - np.sin(theta1)) * (radSsqr / 2.0)

    return area1 - area2
               
def overlap_area(Sx, Sy, Sr, Px, Py, Pr):
    dist = np.hypot(Px - Sx, Py - Sy)

    area = np.piecewise( dist, \
                        [dist <= Sr-Pr, dist >= Sr+Pr], \
                        [np.pi * Pr ** 2, 0, lambda dist: overlap_A(dist, Sr, Pr)])
   
    return area

def dimming(Sx, Sy, Sr, Px, Py, Pr, u):
    dist = np.hypot(Px - Sx, Py - Sy)

    dim = np.piecewise( dist, \
                        [dist >= Sr], \
                        [1, lambda dist: ilim(Sr, dist, u)])
   
    return dim

def overlap_dim( PlanetX, Pr, u, TimeScale):
    Sx = 0
    Sy = 0
    Sr = 1
    Py = 0
    xpos = PlanetX * TimeScale
   
    return ratio_sun_vis( overlap_area( Sx, Sy, Sr, xpos, Py, Pr ) * dimming(Sx, Sy, Sr, xpos, Py, Pr, u ), Sr )


PosPlan = np.linspace(-0.15, 0.15, num=86)

Lum0 = ratio_sun_vis( overlap_area( Sx, Sy, Sr, TimeScale * PosPlan, Py, Pr ), Sr)

#Dim = overlap_dim( PosPlan, Pr, u, TimeScale )

popt, pcov = curve_fit(overlap_dim, Time, relative_flux, p0=[0.125, 0.6, 17])
print(popt)
Dim_opt = overlap_dim( PosPlan, popt[0], popt[1], popt[2] )
perr = np.sqrt(np.diag(pcov))
print('std dev ',': Pr=%g, u=%g, TimeScale=%g' % tuple(perr))

def correlation_func(covariance):
    perr= np.sqrt(np.diag(covariance))
    outer_perr=np.outer(perr,perr)
    correlation = covariance/outer_perr
    correlation[covariance==0]=0
    return correlation

Dim_err_high=overlap_dim( PosPlan, popt[0], popt[1]+perr[1], popt[2] )
Dim_err_low= overlap_dim( PosPlan, popt[0], popt[1]-perr[1], popt[2] )

 

#print('Correlation: ',correlation_func(pcov))
chi_square1 = ((relative_flux - Lum0)**2).sum()
chi_square2 = ((relative_flux - Dim_opt)**2).sum()
print ('Chi Square A ', '= %f' % chi_square1)
print ('Chi Square B ', '= %f' % chi_square2)


plt.figure(1)
#plt.plot(PosPlan, Dim,'go', label='area of circle')
plt.errorbar(Time,relative_flux, yerr=error, fmt= 'bo', label='Flux of star')
plt.plot(PosPlan,Lum0,'ro', label='Area of the circle ')
plt.plot(PosPlan,Dim_opt,'go', label='Limb Darkening')
leg = plt.legend()
plt.xlabel('Days')
plt.ylabel('Relative Flux')
plt.title('Relative Flux vs Days')
plt.grid()

plt.figure(2)
plt.plot(PosPlan,Dim_opt, color='blue', label='optimised u')
plt.plot(PosPlan,Dim_err_high, color='Lime',label='u+1std')
plt.plot(PosPlan,Dim_err_low,color='yellow', label='u-1std')
plt.ylim(.9825,.9900)
plt.xlim(-.075,.075)
leg = plt.legend()
plt.xlabel('Days')
plt.ylabel('Relative Flux')
plt.title('Differant models based of the uncertanty in u')
plt.grid()

Sr=1.1*696.34e6
Pr=1.1*696.34e6*popt[0]
Time_in_sec=(2/18)*24*60*60
TIS_u=perr[2]*24*60*60
Pr_u=perr[0]*Sr
print(Pr,Pr_u)
print('Time of P across dip:',Time_in_sec,'sec +/-',TIS_u)
distance_in_meters=2*Sr
print('Distance of dip:',distance_in_meters,'m')
Pvr=distance_in_meters/Time_in_sec
Pvr_u=(((1/Time_in_sec)**2)*(0)+((distance_in_meters/(Time_in_sec)**2)**2)/((TIS_u)**2))**(1/2)
print('Radial speed of planet',Pvr,'m/sec +/-',Pvr_u)

Svr= 90
Svrmax=110
Svrmin=70
svru=20
mS=1.989e30*1.1
mP=(Svr/Pvr)*mS
mP_u=(((mS/Pvr)**2)*((20)**2)+(((Svr*mS)/(Pvr**2))**2)*((Pvr_u)**2))**(1/2)
print('Mass of planet',mP,'kg +/-',np.abs(mP_u))

G=6.67248e-11
P=((mP+mS)*2*np.pi*G)/((Svr+Pvr)**3)
P_u=(((((2*G*np.pi)/((Svr+Pvr)**3))**2)*(mP_u**2))+((((6*np.pi*G)/((Svr+Pvr)**4))**2)*(svru**2))+((((6*np.pi*G)/((Svr+Pvr)**4))**2)*(Pvr_u**2)))**(1/2)
print('The period is:', P, 's or', P_u,'days +/-',(P_u+P)/(24*60*60))

a=(P*(Pvr+Svr))/(2*np.pi)
a_u=(((((Pvr+Svr)/(np.pi*2))**2)*(P_u**2))+(((P/(np.pi*2))**2)*((Pvr_u)**2))+(((P/(np.pi*2))**2)*((svru)**2)))**(1/2)
print('The combined semi-major axsis is:', a,'m +/-',a_u,'or', a/Sr,'in the Stars Radius +/-', (a_u)/Sr)
Sa=(((mS+mP)*G*(P)**2)/((4*(np.pi)**2)*(1+(mS/mP))**3))**(1/3)
d= (4**(2/3)*(mP**2)*(G**(1/3))*(((mP*3)*(mP+mS)**3)+(mP**4)*(((mS/mP)+1)**3)*(((np.pi**2)*((mS+mP)**2))**(2/3))*(P**2))/(12*(np.pi**2)*((mP**3)**(2/3))*((mP+mS)**6)*((P**2)**(2/3))))
c=(4**(2/3)*(mP**3)*(G**(1/3))*P*(((np.pi**2)*((mS+mP)**2))**(2/3)))/(6*(np.pi**2)*(((mP)**3)**(2/3))*((P**2)**(2/3))*((mP+mS)**2))
Sa_u=((d**2)*(mP_u**2)+((c**2)*(P_u**2)))**(1/2)
Pa=(mS/mP)*Sa
Pa_u=((Sa_u**2)+(a_u**2))**(1/2)
print('star Semi major axsis:',Sa,'m +/-', Sa_u )
print('planet Semi major axsis:',Pa,'m +/-',Pa_u,'or', Pa/Sr,'Radius of the star +/-',Pa_u/Sr)


D_u=(((1/((4/3)*np.pi*(Pr**3)))**2)*(mP_u**2)+(((9*mP)/(4*np.pi*Pr**4))**2)*(Pr_u**2))**(1/2)
print('Density:',mP/((4/3)*np.pi*(popt[0]*Sr)**3),'kg/m^3 +/-',D_u)


print('The arcsin(1/7) is',np.arcsin(1/7.5)*180/np.pi,'Degrees')
py80=(Pa/Sr)*np.cos(80/180*np.pi)
py82=(Pa/Sr)*np.cos(82/180*np.pi)
py84=(Pa/Sr)*np.cos(84/180*np.pi)
py86=(Pa/Sr)*np.cos(86/180*np.pi)
py88=(Pa/Sr)*np.cos(88/180*np.pi)
py90=(Pa/Sr)*np.cos(90/180*np.pi)

def overlap_dim_y( PlanetX, Pr, u, TimeScale, Py):
    Sx = 0
    Sy = 0
    Sr = 1
    xpos = PlanetX * TimeScale
   
    return ratio_sun_vis( overlap_area( Sx, Sy, Sr, xpos, Py, Pr ) * dimming(Sx, Sy, Sr, xpos, Py, Pr, u ), Sr )

Dim_opt_y80 = overlap_dim_y( PosPlan, popt[0], popt[1], popt[2], py80)
Dim_opt_y82 = overlap_dim_y( PosPlan, popt[0], popt[1], popt[2], py82)
Dim_opt_y84 = overlap_dim_y( PosPlan, popt[0], popt[1], popt[2], py84)
Dim_opt_y86 = overlap_dim_y( PosPlan, popt[0], popt[1], popt[2], py86)
Dim_opt_y88 = overlap_dim_y( PosPlan, popt[0], popt[1], popt[2], py88)
Dim_opt_y90 = overlap_dim_y( PosPlan, popt[0], popt[1], popt[2], py90)

plt.figure(3)
plt.plot(PosPlan, Dim_opt_y80, label='i=80')
plt.plot(PosPlan, Dim_opt_y82, label='i=82')
plt.plot(PosPlan, Dim_opt_y84, label='i=84')
plt.plot(PosPlan, Dim_opt_y86, label='i=86')
plt.plot(PosPlan, Dim_opt_y88, label='i=88')
plt.plot(PosPlan, Dim_opt_y90, label='i=90', color='gray')
leg = plt.legend()
plt.xlabel('Days')
plt.ylabel('Relative Flux')
plt.title('Relative Flux based on inclination')
plt.grid()


#plt.figure(3)
#plt.plot(PosPlan,overlap_dim(PosPlan,popt[0],popt[1],18))

