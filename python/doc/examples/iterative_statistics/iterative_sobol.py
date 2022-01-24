#!/usr/bin/env python
# coding: utf-8

# # Computing Sobol' indices with iterative algorithms
# 
# In this example, we show how to estimate Sobol' sensitivity indices using various iterative algorithms. 

# In[1]:


import openturns as ot
from math import sqrt


# We first define the g function. 

# In[2]:


def functionCrue(X) :
    Hd = 3.0
    Zb = 55.5
    L = 5.0e3
    B = 300.0
    Zd = Zb + Hd
    Q, Ks, Zv, Zm = X
    alpha = (Zm - Zv)/L
    H = (Q/(Ks*B*sqrt(alpha)))**(3.0/5.0)
    Zc = H + Zv
    S = Zc - Zd
    return [S]


# In[3]:


g = ot.PythonFunction(4, 1, functionCrue)
g = ot.MemoizeFunction(g)


# Then we create the input distribution. 

# In[4]:


myParam = ot.GumbelMuSigma(1013., 558.)
print(myParam)
Q = ot.ParametrizedDistribution(myParam)
otLOW = ot.TruncatedDistribution.LOWER
Q = ot.TruncatedDistribution(Q, 0, otLOW)
Ks = ot.Normal(30.0, 7.5)
Ks = ot.TruncatedDistribution(Ks, 0, otLOW)
Zv = ot.Uniform(49.0, 51.0)
Zm = ot.Uniform(54.0, 56.0)


# In[5]:


distribution = ot.ComposedDistribution([Q, Ks, Zv, Zm])


# ## Estimating Sobol' indices with a classical algorithm

# We now use the `SobolIndicesExperiment` class in order to create the design of experiments required to estimate Sobol' indices. 

# In[6]:


size = 50
experiment = ot.SobolIndicesExperiment(distribution, size)
inputDesign = experiment.generate()
outputDesign = g(inputDesign)


# Finally, we estimate the Sobol' indices. The first way is to use the classical estimators, which use the full sample. 

# In[7]:


myClassicalSobolStudy = ot.SaltelliSensitivityAlgorithm(inputDesign, outputDesign, size)


# ## The iterative way with Saltelli estimator on a sample

# The second way is to use iterative estimators. In this case, we use the `incrementIndices` method to update the statistics. In order to make accurate computations, we compute the mean of the sample and update the indices with the centered sample instead of the sample itself. 

# In[8]:


myIterativeSobolStudy = ot.SaltelliSobolIndices(g.getInputDimension(), g.getOutputDimension())
muY = outputDesign.computeMean()
centeredOutput = outputDesign - muY
myIterativeSobolStudy.incrementIndices(centeredOutput)


# In order to compare the results, we analyse the first order indices. 

# In[9]:


print(myClassicalSobolStudy.getFirstOrderIndices())
print(myIterativeSobolStudy.getFirstOrderIndices())


# Then we analyse the total order indices. 

# In[10]:


print(myClassicalSobolStudy.getTotalOrderIndices())
print(myIterativeSobolStudy.getTotalOrderIndices())


# ## Iteratively, with Martinez estimator
# 
# The main goal of the iterative estimators is to being able to update the estimators iteratively. We show this with Martinez estimator, which has a `incrementIndices` method that can be used to update the statistics when new points are available. 

# In order to compare the results, we use three algorithms: the first is the classical algorithm, which use the two samples in input. 

# In[11]:


size = 50
experiment = ot.SobolIndicesExperiment(distribution, size)
inputDesign = experiment.generate()
outputDesign = g(inputDesign)


# In[12]:


myClassicalSobolStudy = ot.MartinezSensitivityAlgorithm(inputDesign, outputDesign, size)


# Then we create an iterative algorithm and update the statistics with sub-samples of smaller size. Notice that the `inputDesign` sample in the nex `for` loop has a size equal to $p+2 = 6$, where $p$ is the input dimension. 

# In[13]:


myIterativeSobolStudy = ot.MartinezSobolIndices(g.getInputDimension(), g.getOutputDimension())


# In[27]:


for i in range(size):
    experiment = ot.SobolIndicesExperiment(distribution, 1)
    inputDesign = experiment.generate()
    outputDesign = g(inputDesign)
    myIterativeSobolStudy.incrementIndices(outputDesign)


# The last way is to use the all sample in one row. This will perform iterations from the sample. 

# In[15]:


size = 50
experiment = ot.SobolIndicesExperiment(distribution, size)
inputDesign = experiment.generate()
outputDesign = g(inputDesign)


# In[16]:


myIterativeSobolStudy2 = ot.MartinezSobolIndices(g.getInputDimension(), g.getOutputDimension())


# In[17]:


myIterativeSobolStudy2.incrementIndices(outputDesign)


# In[18]:


print("First order")
print(myClassicalSobolStudy.getFirstOrderIndices())
print(myIterativeSobolStudy.getFirstOrderIndices())
print(myIterativeSobolStudy2.getFirstOrderIndices())
print("Total order")
print(myClassicalSobolStudy.getTotalOrderIndices())
print(myIterativeSobolStudy.getTotalOrderIndices())
print(myIterativeSobolStudy2.getTotalOrderIndices())


# ## The iterative way with Jansen estimator

# In[19]:


size = 50
experiment = ot.SobolIndicesExperiment(distribution, size)
inputDesign = experiment.generate()
outputDesign = g(inputDesign)


# In[20]:


myClassicalSobolStudy = ot.JansenSensitivityAlgorithm(inputDesign, outputDesign, size)
myIterativeSobolStudy = ot.JansenSobolIndices(g.getInputDimension(), g.getOutputDimension())
myIterativeSobolStudy.incrementIndices(outputDesign)


# In[21]:


print("First order")
print(myClassicalSobolStudy.getFirstOrderIndices())
print(myIterativeSobolStudy.getFirstOrderIndices())
print("Total order")
print(myClassicalSobolStudy.getTotalOrderIndices())
print(myIterativeSobolStudy.getTotalOrderIndices())


# ## Iteratively with Mauntz-Kucherenko algorithm

# In[22]:


size = 50
experiment = ot.SobolIndicesExperiment(distribution, size)
inputDesign = experiment.generate()
outputDesign = g(inputDesign)


# In[23]:


myClassicalSobolStudy = ot.MauntzKucherenkoSensitivityAlgorithm(inputDesign, outputDesign, size)


# In[24]:


myIterativeSobolStudy = ot.MauntzKucherenkoSobolIndices(g.getInputDimension(), g.getOutputDimension())
muY = outputDesign.computeMean()
centeredOutput = outputDesign - muY
myIterativeSobolStudy.incrementIndices(centeredOutput)


# In[25]:


print("First order")
print(myClassicalSobolStudy.getFirstOrderIndices())
print(myIterativeSobolStudy.getFirstOrderIndices())
print("Total order")
print(myClassicalSobolStudy.getTotalOrderIndices())
print(myIterativeSobolStudy.getTotalOrderIndices())

