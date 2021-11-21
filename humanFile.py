#!/usr/bin/env python
# coding: utf-8

# In[5]:


from pyopenms import *
from urllib.request import urlretrieve
dig = ProteaseDigestion()
dig.getEnzymeName() 
hum = "".join([l.strip() for l in open("uniprot-yourlist_M202111204ABAA9BC7178C81CEBC9459510EDDEA3300C607.fasta").readlines()[1:]])
hum = AASequence.fromString(hum)

result = []
dig.digest(hum, result)
print(result[4].toString())
len(result) 


# In[ ]:




