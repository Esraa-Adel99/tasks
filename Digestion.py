#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyopenms import *
from urllib.request import urlretrieve
gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
urlretrieve (gh + "/src/data/P02769.fasta", "bsa.fasta")

dig = ProteaseDigestion()
dig.getEnzymeName() 
bsa = "".join([l.strip() for l in open("bsa.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
result = []
dig.digest(bsa, result)
print(result[4].toString())
len(result)


# In[3]:


dig.digest(bsa, result, 7, 40)
for s in result:
    print(s.toString())


# In[4]:


dig.setMissedCleavages(2)
dig.digest(bsa, result, 7, 40)
for s in result:
    print(s.toString())


# In[5]:


names = []
ProteaseDB().getAllNames(names)
len(names) 
e = ProteaseDB().getEnzyme('Lys-C')
e.getRegExDescription()
e.getRegEx()


# In[6]:


from urllib.request import urlretrieve
gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
urlretrieve (gh + "/src/data/P02769.fasta", "bsa.fasta")
dig = ProteaseDigestion()
dig.setEnzyme('Lys-C')
bsa = "".join([l.strip() for l in open("bsa.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
result = []
dig.digest(bsa, result)
print(result[4].toString())
len(result)


# In[ ]:




