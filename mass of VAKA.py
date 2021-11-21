#!/usr/bin/env python
# coding: utf-8

# In[14]:


from pyopenms import *
seq = AASequence.fromString("VAKA")
mfull = seq.getMonoWeight()
print("Monoisotopic mass of peptide [M] is", mfull)
print("*********************************")
totalmass = 0
print("The peptide", str(seq), "consists of the following amino acids:")
for aa in seq:
    print(aa.getName(), ":", aa.getMonoWeight())
    totalmass +=aa.getMonoWeight()
    
print("*********************************")    
print("mass(V) + mass(A) + mass(K) + mass(A) = ",totalmass)    


# In[ ]:




