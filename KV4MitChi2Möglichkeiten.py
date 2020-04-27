#!/usr/bin/env python
# coding: utf-8

# In[3]:


import scipy
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.stats import linregress
import numpy as np
import scipy.fftpack
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab
import math as m



#Gain = mu1-mu0 = mupeak - mupedestal
#Std = varianz**2 (in fit ist standardabweichung)
def trennlinie():
    print('______________________________________________________________________________')
def gaussian(x, a, mu, sigma):
    return((a/(np.sqrt(2.*np.pi)*sigma**2))*np.exp((-1/2)*((x-mu)/sigma)**2))
def poisson(f0, mu):
    return(((mu**f0)*m.exp(-mu))/m.factorial(f0))
def F():
    return N/Ntot
#def mu(G,P0,i):
#    return P0 + i*G
#def mu(data):
#    return np.average(data)

#def sigma(data):
#   return np.std(data)
#________________________________________________________________________
print("4.3.2 & 4.3.3")
#LED 54 Volt
#Für Pedelstal
print("LED 64 Volt")
data = np.loadtxt("C:/Users/jloem/Desktop/KurzversuchAblage/KV4 Messungen/LED.64V.txt", skiprows = 8)
pedestalX = data[:,0][(data[:,0]>-27)&(data[:,0]<27)]
pedestalY = data[:,1][(data[:,0]>-27)&(data[:,0]<27)]
var, cov = curve_fit(gaussian, pedestalX, pedestalY)
plt.plot(pedestalX,pedestalY)
plt.plot(pedestalX, gaussian(pedestalX, var[0], var[1], var[2]))
pylab.xlim([-50,100])
plt.yscale('log')
plt.show()

#Zahl der Ereignisse: f0
f0 = (np.sum(pedestalY))
mu = (var[1])
print("Mittelwert(Pedestal)", var[1], "Standardabweichung:", var[2])
print("Zahl der Ereignisse im Pedestel:", f0)

#Für ersten Peak
peakX = data[:,0][(data[:,0]>272)&(data[:,0]<350)]
peakY = data[:,1][(data[:,0]>272)&(data[:,0]<350)]
varp, covp = curve_fit(gaussian, peakX, peakY, p0=(15000, 310, 10))
plt.plot(peakX,peakY )
plt.plot(peakX, gaussian(peakX, varp[0], varp[1], varp[2]))
pylab.xlim([260,360])
plt.yscale('log')
plt.show()


gain = varp[1]-var[1]
print("Mittelwert(Erseter Peak):", varp[1], "Standardabweichung:", varp[2])
print("Gain:", gain)
#print("als relativer Fehler auf den Mittelwert ergbit sich:", (1/m.sqrt(varp[1])) KP warum das nicht klappt damit...


#__________________________________________________________________


#Unter Annahme einer Poissonverteilung warsch F() = N0.5/Ntot = exp(-mu) mu = varianz = (std)**2
#64 Volt f>0.5
dataX = data[:,0][(data[:,0]>(gain)/2)&(data[:,0]<1000)]
dataY = data[:,1][(data[:,0]>(gain)/2)&(data[:,0]<1000)]
N = np.sum(dataY)

#64 Volt ftot
DataX = data[:,0][(data[:,0]>-27)&(data[:,0]<1000)]
DataY = data[:,1][(data[:,0]>-27)&(data[:,0]<1000)]
Ntot = np.sum(DataY)

Mittelwerp = -np.log(F())
relFehler = (1/m.sqrt(Mittelwerp))

print("Unter Annahme einer Poissonverteilung ist der Mittelwert gleich der Varianz:", Mittelwerp)
print("als relativer Fehler auf den Mittelwert ergbit sich:", relFehler)
#_______________________________________________________________
#Für Histogramm
histoX = data[:,0][(data[:,0]>-27)&(data[:,0]<1000)]
histoY = data[:,1][(data[:,0]>-27)&(data[:,0]<1000)]
plt.hist(histoY, bins = 'auto', range=(-10,1000))#könnte noch normieren bzw bins anpassen
plt.title("Histogramm LED 64 Volt")
plt.xlabel("Counts")
plt.ylabel("Events")
plt.yscale('log')
plt.show()

Mittelwert = sum((histoY*histoX)/Ntot)/gain
print("Aus dem Histogramm, Mittelwert:", Mittelwert, "Standardabweichung:", histoY.std())# für alle Histogramme berechnen!


# In[32]:


print("4.2.3 & 4.2.4")
#________________________________________________________

##Events oberhalt f1/2
def F():
    return N/Ntot
#DCR
def DCR():
    return F()/(144*10**(-9))
#CrossTalkWarscheinlichkeit n für f>1.5
def Pcn():
    return n/Ntot
#________________________________________________________
#Dunkelspektren 4.2.3
#DCR 56 Volt   
print("Dark Count Rate 56 Volt")
dark = np.loadtxt("C:/Users/jloem/Desktop/KurzversuchAblage/KV4 Messungen/DCR.56.03V_histo.txt", skiprows = 8)
#Gatelänge 144ns
print("Gatelänge 144ns")
#Gain
print("Gain:", gain)
#DCR 56 Volt N>0.5
darkX5 = dark[:,0][(dark[:,0]>56)&(dark[:,0]<1000)]
darkY5 = dark[:,1][(dark[:,0]>56)&(dark[:,0]<1000)]
N = np.sum(darkY5)

#DCR 56 Volt n>1.5
darkPX = dark[:,0][(dark[:,0]>168)&(dark[:,0]<1000)]
darkPY = dark[:,1][(dark[:,0]>168)&(dark[:,0]<1000)]
n = np.sum(darkPY)

#DCR 56 Volt gesamt
darkX = dark[:,0][(dark[:,0]>-56)&(dark[:,0]<1000)]
darkY = dark[:,1][(dark[:,0]>-56)&(dark[:,0]<1000)]
Ntot = np.sum(darkY)

plt.plot(darkX,darkY)
plt.yscale('log')

print("DCR:", DCR(), "[Hz]")
print("Cross-Talk Warscheinlichkeit", Pcn())
#________________________________________________________
trennlinie()


# In[34]:


#Gain und Spannung
print("4.3.4")
x=[54, 56, 58, 60, 62, 64]
X = sum(x)-np.mean(x)
#[gain1, gain2, gain3, gain4, gain5, gain6] #Werte Eintragen
y=[59.01613105267129, 112.55272551330273, 162.17202375117805, 214.74143168392732, 261.03739357038853, 312.22293918485394]
pylab.xlim([50,70])
plt.xlabel("Spannung [Volt]")
plt.ylabel("Gain")
plt.plot(x,y)
plt.show()
#_____________________________________________________________

#Linearer Fit
z = np.polyfit(x, y, 1)
f = np.poly1d(z)
plt.plot((0, max(x)), (f(0), f(max(x))), 'r')

#_____________________________________________________________
# Das Achsen-Objekt des Diagramms in einer Variablen ablegen:
ax = plt.gca()

# Die obere und rechte Achse unsichtbar machen:
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

# Die linke Diagrammachse auf den Bezugspunkt '0' der x-Achse legen:
ax.spines['left'].set_position(('data',0))

# Die untere Diagrammachse auf den Bezugspunkt '0' der y-Achse legen:
ax.spines['bottom'].set_position(('data',0))

# Ausrichtung der Achsen-Beschriftung festlegen:
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
#_____________________________________________________________

plt.ylabel("Gain")
plt.xlabel("Spannung [Volt]")
print("Extrapolation")
plt.show()
#Für Vergleich der Werte
print("Fit(x)=",f)

#Fehler und StandardAbweichung
linregress(x,y)


# In[4]:


#Model___________________________________!!!!!! alle peaks werden nach Gauss gefittet, die höhen dann nach poisson!!!!!!
print("Modell SiPM")
import scipy
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.stats import linregress
import numpy as np
import scipy.fftpack
import scipy.signal
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pylab
import math as m

def fact(x):
	for i in x:
		return m.gamma(i + 1.0)

def poisson(i, mu, h):
	return ((mu**i)*m.exp(-mu))/fact(i)*h

def normal(x, mu, sigma):
	return 1.0/(np.sqrt(2.0*np.pi*sigma*sigma))*np.exp(-0.5*((x-mu)/sigma)**2)

def fitfunction(x, mu, sigma, h):
	return h*normal(x, mu, sigma)

def plotSpectra(data):
    x = data[:,0][(data[:,0]>-27)&(data[:,0]<1000)] #____________XY um 100% Array zu haben
    y = data[:,1][(data[:,0]>-27)&(data[:,0]<1000)]
    plt.plot(x,y)
    plt.show()

def readData(_filename):
    print('Read data from',_filename)
    f = open(_filename)
    data = np.loadtxt(f, skiprows = 8) #___________array
    return data

def fitSpectrum(x, y, startMu):
    var, cov = curve_fit(fitfunction, x, y, p0 =(startMu, 100.0, 1.0))
    plt.plot(x, fitfunction(x, var[0], var[1], var[2]))
    #plt.yscale('log')
    #plt.show()
    print("mu=%f ± %f, sigma=%f ± %f, h=%f ± %f" % (var[0], cov[0,0], var[1], cov[1,1], var[2], cov[2,2]))
    return (var[0], var[1], var[2])

def normalizeHistogram(y):
	s = np.sum(y)
	return y/s

def main(_filename):
    data = readData(_filename)
    plotSpectra(data)
    x = data[:,0][(data[:,0]>-27)&(data[:,0]<1000)]
    y = data[:,1][(data[:,0]>-27)&(data[:,0]<1000)] 
    y = normalizeHistogram(y)
    plt.plot(x,y)
    Is = np.array([0, 1, 2, 3])
    Xs = np.array([0.1, 1.1, 2.1, 3.1]) #Ist das array, dass das Intervall zum jeweiligen fit gibt
                                            #Wirkürliche Kommazahlen um das Array zu initialisieren
    Sigmas = np.array([0.1, 1.1, 2.1, 3.1]) #array mit jeweiligem startfehler
    Hs = np.array([0.1, 1.1, 2.1, 3.1])# 
    print("Peak 0:")
    intv = x < 100	# Select intervall for peak
    Xs[0], Sigmas[0], Hs[0] = fitSpectrum(x[intv],y[intv],0.0)#jeweils um Gain erweitert
    print("Peak 1:")
    intv = (x>200)&(x<450) 	# Select intervall for peak
    Xs[1], Sigmas[1], Hs[1] = fitSpectrum(x[intv],y[intv],300.0)
    print("Peak 2:")
    intv = (x>500)&(x<750) 	# Select intervall for peak
    Xs[2], Sigmas[2], Hs[2] = fitSpectrum(x[intv],y[intv],620.0)
    print("Peak 3:")
    intv = (x > 800)&(x<1050) 	# Select intervall for peak
    Xs[3], Sigmas[3], Hs[3] = fitSpectrum(x[intv],y[intv],920.0)
    #___________________________________________________________________Chi**2 test
    print("Methode der kleinsten Quadrate")
    
    chiqua1 = np.sum(((Xs[0] -0.592561)/Sigmas[0])**2)
    chiqua2 = np.sum(((Xs[1] -312.817051)/Sigmas[1])**2)
    chiqua3 = np.sum(((Xs[2] -625.185958)/Sigmas[2])**2)
    chiqua4 = np.sum(((Xs[3] -936.454221)/Sigmas[3])**2)
    print("Mit Sigma:",chiqua1, chiqua2, chiqua3, chiqua4)
    chiqua1 = np.sum(((Xs[0] -0.592561)/0.000197)**2)
    chiqua2 = np.sum(((Xs[1] -0.592561)/0.006509)**2)
    chiqua3 = np.sum(((Xs[2] -625.185958)/0.073602)**2)
    chiqua4 = np.sum(((Xs[3] -936.454221)/0.305927)**2)
    print("Mit Fehler auf Sigma:",chiqua1, chiqua2, chiqua3, chiqua4)
    chiqua1 = np.sum(((Xs[0] -0.592561)/0.592561)**2)
    chiqua2 = np.sum(((Xs[1] -312.817051)/312.817051)**2)
    chiqua3 = np.sum(((Xs[2] -625.185958)/625.185958)**2)
    chiqua4 = np.sum(((Xs[3] -936.454221)/936.454221)**2)
    print("Mit Mittelwert:",chiqua1, chiqua2, chiqua3, chiqua4)
    #____________________________________________________________________
    print("Possion Fit:")
    var, cov = curve_fit(poisson, xdata=Is, ydata=Hs, p0=(0.1, 0.03))
    print("lambda= %f ± %f, Höhe=%f ± %f" % (var[0],cov[0,0],var[1],cov[1,1]))
    plt.plot(Xs, normal(Xs, Xs, Sigmas)*poisson(Is, var[0], var[1]))#die schlechte Kurve
    plt.hist(x, bins=1000, weights=y)#y ist auf 1 normiert(geteilt duch sum(y)), anzahl bins = sum
    plt.show()
   

if __name__=="__main__":
    main("C:/Users/jloem/Desktop/KurzversuchAblage/KV4 Messungen/LED.64V.txt")

