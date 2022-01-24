#!/usr/bin/env python
# coding: utf-8

# # Computing the mean or variance with an iterative algorithm
# 
# In this example, we show how to use the `IterativeMean` and `IterativeVariance` classes. These algorithms compute an estimate of the mean and variance using iterative algorithms. In other words, these algorithm can be used to iteratively update the estimate of the mean or variance using one point at a time instead of using the full sample. 

# In[1]:


import openturns as ot


# We first create the two estimators with given dimension. 

# In[2]:


dim = 5
myMean = ot.IterativeMean(dim)
myVariance = ot.IterativeVariance(dim)


# Then we run the simulations. In this example, we generate 50 points from a gaussian distriution. For each point, we update the mean and variance using the `increment` method. 

# In[3]:


n = ot.Normal(dim)
size = 50
# 2.1. Increment with several points
for i in range(size):
    point = n.getRealization()
    myMean.increment(point)
    myVariance.increment(point)


# We can now see the results of the algorithm. The `getIteration` method returns the sample size used in the estimation. Since estimating the variance requires to estimate the mean, the variance has a `getMean` method too. Hence `myMean.getMean()` and `myVariance.getMean()` should be the same. 

# In[4]:


print ("Iteration " + str(myVariance.getIteration()))
print (myMean.getMean())
print (myVariance.getMean())
print (myVariance.getStandardDeviation())
print (myVariance.getVariance())


# Although the main point of the iterative algorithm is to update the statistics with one single point, we can also update the statistics with a `Sample`. This might be useful in the case where several points are available at the same time, e.g. with parallel computing.  

# In[5]:


sample = n.getSample(size)
myMean.increment(sample)
myVariance.increment(sample)


# Notice that two sets of 50 points have been used, leading to a total of 100 points. 

# In[6]:


print ("Iteration " + str(myVariance.getIteration()))
print (myMean.getMean())
print (myVariance.getMean())
print (myVariance.getStandardDeviation())
print (myVariance.getVariance())

