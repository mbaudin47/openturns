#!/usr/bin/env python
# coding: utf-8

# # Iterative mean
# 
# In this example, we show how to use the `IterativeMean` class. This algorithm computes an estimate of the mean using iterative algorithms. In other words, this algorithm can be used to iteratively update the estimate of the mean using one point at a time instead of using the full sample. 

# In[1]:


import openturns as ot


# Create the mean object.

# In[2]:


dim = 5
myMean = ot.IterativeMean(dim)


# Perform the simulations. In the following session, we increment using the `increment` method one `Point` at a time.

# In[3]:


n = ot.Normal(dim)
size = 50
for i in range(size):
    point = n.getRealization()
    myMean.increment(point)


# The `getMean` method returns the sample mean.

# In[4]:


myMean.getMean()


# The `getIteration` method return the sample size.

# In[5]:


myMean.getIteration()


# We can also increment with a `Sample`.

# In[6]:


sample = n.getSample(size)
myMean.increment(sample)


# In[7]:


print ("Iteration " + str(myMean.getIteration()))
print (myMean.getMean())

