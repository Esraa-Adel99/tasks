#!/usr/bin/env python
# coding: utf-8

# In[3]:


from pyopenms import *
from urllib.request import urlretrieve
dig = ProteaseDigestion()
dig.getEnzymeName() 
bsa = "".join([l.strip() for l in open("uniprot-yourlist_M20211120A084FC58F6BBA219896F365D15F2EB442FD6B5M.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)

result = []
dig.digest(bsa, result)
print(result[4].toString())
len(result) 


# In[ ]:




