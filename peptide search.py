#!/usr/bin/env python
# coding: utf-8

# In[1]:


from urllib.request import urlretrieve
from pyopenms import *
gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
urlretrieve (gh +"/src/data/SimpleSearchEngine_1.mzML", "searchfile.mzML")
urlretrieve (gh +"/src/data/SimpleSearchEngine_1.fasta", "search.fasta")
protein_ids = []
peptide_ids = []
SimpleSearchEngineAlgorithm().search("searchfile.mzML", "search.fasta", protein_ids, peptide_ids)


# In[2]:


for peptide_id in peptide_ids:
 
  print (35*"=")
  print ("Peptide ID m/z:", peptide_id.getMZ())
  print ("Peptide ID rt:", peptide_id.getRT())
  print ("Peptide scan index:", peptide_id.getMetaValue("scan_index"))
  print ("Peptide scan name:", peptide_id.getMetaValue("scan_index"))
  print ("Peptide ID score type:", peptide_id.getScoreType())
  
  for hit in peptide_id.getHits():
    print(" - Peptide hit rank:", hit.getRank())
    print(" - Peptide hit charge:", hit.getCharge())
    print(" - Peptide hit sequence:", hit.getSequence())
    mz = hit.getSequence().getMonoWeight(Residue.ResidueType.Full, hit.getCharge()) / hit.getCharge()
    print(" - Peptide hit monoisotopic m/z:", mz)
    print(" - Peptide ppm error:", abs(mz - peptide_id.getMZ())/mz *10**6 )
    print(" - Peptide hit score:", hit.getScore())


# In[3]:


tsg = TheoreticalSpectrumGenerator()
thspec = MSSpectrum()
p = Param()
p.setValue("add_metainfo", "true")
tsg.setParameters(p)
peptide = AASequence.fromString("RPGADSDIGGFGGLFDLAQAGFR")
tsg.getSpectrum(thspec, peptide, 1, 1)

for ion, peak in zip(thspec.getStringDataArrays()[0], thspec):
    print(ion, peak.getMZ())

e = MSExperiment()
MzMLFile().load("searchfile.mzML", e)
print ("Spectrum native id", e[2].getNativeID() )
mz,i = e[2].get_peaks()
peaks = [(mz,i) for mz,i in zip(mz,i) if i > 1500 and mz > 300]
for peak in peaks:
  print (peak[0], "mz", peak[1], "int")


# In[4]:


salgo = SimpleSearchEngineAlgorithm()
p = salgo.getDefaults()
print ( p.items() )
p[b'precursor:mass_tolerance'] = 4.0
salgo.setParameters(p)

protein_ids = []
peptide_ids = []
salgo.search("searchfile.mzML", "search.fasta", protein_ids, peptide_ids)
print("Found", len(peptide_ids), "peptides")


# In[ ]:





# In[ ]:




