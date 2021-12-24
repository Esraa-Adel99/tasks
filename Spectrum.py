#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyopenms import *
spectrum = MSSpectrum()
mz = range(1500, 500, -100)
i = [0 for mass in mz]
spectrum.set_peaks([mz, i])
spectrum.sortByPosition()
for p in spectrum:
  print(p.getMZ(), p.getIntensity())
for mz, i in zip(*spectrum.get_peaks()):
  print(mz, i)
print(spectrum[2].getMZ(), spectrum[2].getIntensity())


# In[2]:


spectrum = MSSpectrum()
spectrum.setDriftTime(25) 
spectrum.setRT(205.2) 
spectrum.setMSLevel(3) 
spectrum.set_peaks( ([401.5], [900]) )
p = Precursor()
p.setMZ(600) 
p.setIsolationWindowLowerOffset(1.5)
p.setIsolationWindowUpperOffset(1.5)
p.setActivationEnergy(40) 
p.setCharge(4) 
spectrum.setPrecursors( [p] )
IS = InstrumentSettings()
IS.setPolarity(IonSource.Polarity.POSITIVE)
spectrum.setInstrumentSettings(IS)
polarity = spectrum.getInstrumentSettings().getPolarity()
if (polarity == IonSource.Polarity.POSITIVE):
  print("scan polarity: positive")
elif (polarity == IonSource.Polarity.NEGATIVE):
  print("scan polarity: negative")
fda = FloatDataArray()
fda.setName("Signal to Noise Array")
fda.push_back(15)
sda = StringDataArray()
sda.setName("Peak annotation")
sda.push_back("y15++")
spectrum.setFloatDataArrays( [fda] )
spectrum.setStringDataArrays( [sda] )
exp = MSExperiment()
exp.addSpectrum(spectrum)
spectrum2 = MSSpectrum()
spectrum2.set_peaks( ([1, 2], [1, 2]) )
exp.addSpectrum(spectrum2)
MzMLFile().store("testfile.mzML", exp)


# In[6]:


import matplotlib.pyplot as plt

def plot_spectrum(spectrum):
    
    for mz, i in zip(*spectrum.get_peaks()):
        plt.plot([mz, mz], [0, i], color = 'black')
        plt.text(mz, i, str(mz))

    
    title = ''
    if spectrum.getRT() >= 0:
        title += 'RT: ' + str(spectrum.getRT())
    if len(spectrum.getPrecursors()) >= 1:
        title += '   Precursor m/z: ' + str(spectrum.getPrecursors()[0].getMZ())

    plt.title(title)
    plt.ylabel('intensity')
    plt.xlabel('m/z')
    plt.ylim(bottom=0)

    plt.show()


plot_spectrum(spectrum)


# In[ ]:




