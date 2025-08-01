#!/usr/bin/env python

import openturns as ot
import openturns.testing as ott

ot.TESTPREAMBLE()

# QQPlot tests
size = 100
normal = ot.Normal(1)
sample = normal.getSample(size)
sample2 = ot.Gamma(3.0, 4.0, 0.0).getSample(size)
twoSamplesQQPlot = ot.VisualTest.DrawQQplot(sample, sample2)
print("twoSamplesQQPlot = ", twoSamplesQQPlot)

sampleDistributionQQPlot = ot.VisualTest.DrawQQplot(sample, normal)
print("sampleDistributionQQPlot = ", sampleDistributionQQPlot)

dist = ot.Geometric()
qq_plot = ot.VisualTest.DrawQQplot(dist.getSample(size), dist)
print("discrete QQPlot = ", qq_plot)

sample2 = normal.getSample(size)
twoSamplesPPPlot = ot.VisualTest.DrawPPplot(sample, sample2)
print("twoSamplesPPPlot = ", twoSamplesPPPlot)

sampleDistributionPPplot = ot.VisualTest.DrawPPplot(sample, normal)
print("sampleDistributionPPplot = ", sampleDistributionPPplot)

# HenryLine test
size = 100
normal = ot.Normal(1)
sample = normal.getSample(size)
henryPlot = ot.VisualTest.DrawHenryLine(sample)
print("HenryPlot = ", henryPlot)

# LinearModel tests
dimension = 2
R = ot.CorrelationMatrix(dimension)
R[0, 1] = 0.8
distribution = ot.Normal(ot.Point(dimension, 3.0), ot.Point(dimension, 2.0), R)
size = 100
sample2D = distribution.getSample(size)
firstSample = ot.Sample(size, 1)
secondSample = ot.Sample(size, 1)
for i in range(size):
    firstSample[i] = ot.Point(1, sample2D[i, 0])
    secondSample[i] = ot.Point(1, sample2D[i, 1])

lmtest = ot.LinearModelAlgorithm(firstSample, secondSample).getResult()
drawLinearModelVTest = ot.VisualTest.DrawLinearModel(lmtest)
print("LinearModelV = ", drawLinearModelVTest)

drawLinearModelResidualTest = ot.VisualTest.DrawLinearModelResidual(lmtest)
print("LinearModelR = ", drawLinearModelResidualTest)

# Parallel coordinates tests
size = 100
inputDimension = 6
inputSample = ot.Normal(inputDimension).getSample(size)
inputVar = ["X" + str(i) for i in range(inputDimension)]
formula = ot.Description(1)
expression = ""
for i in range(inputDimension):
    if i > 0:
        expression += "+"
    expression += "cos(" + str(i + 1) + "*" + inputVar[i] + ")"
formula[0] = expression
model = ot.SymbolicFunction(inputVar, formula)
outputSample = model(inputSample)
parallelCoordinatesValue = ot.VisualTest.DrawParallelCoordinates(
    inputSample, outputSample, 2.5, 3.0, "red", False
)
print("parallelCoordinatesValue = ", parallelCoordinatesValue)

parallelCoordinatesQuantile = ot.VisualTest.DrawParallelCoordinates(
    inputSample, outputSample, 0.7, 0.9, "red", False
)
print("parallelCoordinatesQuantile = ", parallelCoordinatesQuantile)

# KendallPlot tests
size = 100
copula1 = ot.FrankCopula(1.5)
copula2 = ot.GumbelCopula(4.5)
sample1 = copula1.getSample(size)
sample1.setName("data 1")
sample2 = copula2.getSample(size)
sample2.setName("data 2")
kendallPlot1 = ot.VisualTest.DrawKendallPlot(sample1, copula2)
print("KendallPlot1 = ", kendallPlot1)

kendallPlot2 = ot.VisualTest.DrawKendallPlot(sample2, sample1)
print("KendallPlot2 = ", kendallPlot2)

# Clouds
sample = ot.Normal(4).getSample(200)
clouds = ot.VisualTest.DrawPairs(sample)
print("Clouds = ", clouds)
distribution = ot.JointDistribution(
    [ot.HistogramFactory().build(sample.getMarginal(i)) for i in range(4)]
)
cloudsMarginals = ot.VisualTest.DrawPairsMarginals(sample, distribution)
print("CloudsMarginals = ", cloudsMarginals)

# dependence functions
copula = ot.GumbelCopula()
data = copula.getSample(100000)
graph1 = ot.VisualTest.DrawUpperTailDependenceFunction(data)
graph2 = ot.VisualTest.DrawUpperExtremalDependenceFunction(data)
graph3 = ot.VisualTest.DrawLowerTailDependenceFunction(data)
graph4 = ot.VisualTest.DrawLowerExtremalDependenceFunction(data)

# check vs theoretical value
theta = copula.getTheta()
ref = 2.0 - 2.0 ** (1.0 / theta)
value = graph1.getDrawable(0).getData()[-4, 1]
ott.assert_almost_equal(value, ref, 1e-2, 1e-3)

# Draw points inside and outside a domain
domain = ot.Domain(ot.Interval([-2.0, 0.0, -1], [2.0, 3.0, 1.0]))
U = ot.Uniform(-3, 3)
dist = ot.JointDistribution([U, U, U])
grid = ot.VisualTest.DrawInsideOutside(domain, dist.getSample(30))
print(grid)

# CDF plot with unique values
dist = ot.Mixture([ot.Uniform(8, 15), ot.Dirac(10)], [0.1, 0.9])
ud = ot.UserDefined(dist.getSample(1000))
graph = ot.VisualTest.DrawCDFplot(ud.getSample(100), ud)
