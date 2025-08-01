#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import math as m


def cleanPoint(inPoint):
    dim = inPoint.getDimension()
    for i in range(dim):
        if m.fabs(inPoint[i]) < 1.0e-10:
            inPoint[i] = 0.0
    return inPoint


def cleanCovariance(inCovariance):
    dim = inCovariance.getDimension()
    for j in range(dim):
        for i in range(dim):
            if m.fabs(inCovariance[i, j]) < 1.0e-10:
                inCovariance[i, j] = 0.0
    return inCovariance


# Instantiate one distribution object
referenceDistribution = [
    ot.TruncatedNormal(2.0, 1.5, 1.0, 4.0),
    ot.TruncatedNormal(2.0, 1.5, 1.0, 200.0),
    ot.TruncatedNormal(2.0, 1.5, -200.0, 4.0),
    ot.TruncatedNormal(2.0, 1.5, 1.0, 4.0),
]
distribution = [
    ot.TruncatedDistribution(ot.Normal(2.0, 1.5), 1.0, 4.0),
    ot.TruncatedDistribution(ot.Normal(2.0, 1.5), 1.0, ot.TruncatedDistribution.LOWER),
    ot.TruncatedDistribution(ot.Normal(2.0, 1.5), 4.0, ot.TruncatedDistribution.UPPER),
    ot.TruncatedDistribution(
        ot.Normal(2.0, 1.5), ot.Interval([1.0], [4.0], [True], [True])
    ),
]

# add a 2-d test
dimension = 2
# This distribution takes too much time for the test
# size = 70
# ref = ot.Normal(dimension)
# sample = ref.getSample(size)
# ks = ot.KernelSmoothing().build(sample)
# Use a multivariate Normal distribution instead
ks = ot.Normal(2)
truncatedKS = ot.TruncatedDistribution(
    ks, ot.Interval([-0.5] * dimension, [2.0] * dimension)
)
distribution.append(truncatedKS)
referenceDistribution.append(ks)  # N/A
# Add a non-truncated example
weibull = ot.WeibullMin(2.0, 3.0)
distribution.append(ot.TruncatedDistribution(weibull))
referenceDistribution.append(weibull)
ot.RandomGenerator.SetSeed(0)

for testCase in range(len(distribution)):
    print("Distribution ", distribution[testCase])

    # Is this distribution elliptical ?
    print("Elliptical = ", distribution[testCase].isElliptical())

    # Is this distribution continuous ?
    print("Continuous = ", distribution[testCase].isContinuous())

    # Test for realization of distribution
    oneRealization = distribution[testCase].getRealization()
    print("oneRealization=", repr(oneRealization))

    # Test for sampling
    size = 10000
    oneSample = distribution[testCase].getSample(size)
    print("oneSample first=", repr(oneSample[0]), " last=", repr(oneSample[size - 1]))
    print("mean=", repr(oneSample.computeMean()))
    print("covariance=", repr(oneSample.computeCovariance()))

    # Define a point
    point = ot.Point(distribution[testCase].getDimension(), 1.5)
    print("Point= ", repr(point))

    # Show PDF and CDF of point
    eps = 1e-5

    DDF = distribution[testCase].computeDDF(point)
    print("ddf      =", repr(DDF))
    print("ddf (ref)=", repr(referenceDistribution[testCase].computeDDF(point)))

    PDF = distribution[testCase].computePDF(point)
    print("pdf      =%.6f" % PDF)
    print("pdf (ref)=%.6f" % referenceDistribution[testCase].computePDF(point))

    CDF = distribution[testCase].computeCDF(point)
    print("cdf=%.6f" % CDF)
    CCDF = distribution[testCase].computeComplementaryCDF(point)
    print("ccdf=%.6f" % CCDF)
    print("cdf (ref)=%.6f" % referenceDistribution[testCase].computeCDF(point))

    PDFgr = distribution[testCase].computePDFGradient(point)
    print("pdf gradient      =", repr(cleanPoint(PDFgr)))

    CDFgr = distribution[testCase].computeCDFGradient(point)
    print("cdf gradient      =", repr(cleanPoint(CDFgr)))

    # quantile
    quantile = distribution[testCase].computeQuantile(0.95)
    print("quantile=", repr(quantile))
    print("quantile=", repr(referenceDistribution[testCase].computeQuantile(0.95)))
    print("cdf(quantile)=%.6f" % distribution[testCase].computeCDF(quantile))
    # Get 95% survival function
    inverseSurvival = ot.Point(
        distribution[testCase].computeInverseSurvivalFunction(0.95)
    )
    print("InverseSurvival=", repr(inverseSurvival))
    print(
        "Survival(inverseSurvival)=%.6f"
        % distribution[testCase].computeSurvivalFunction(inverseSurvival)
    )
    print("entropy=%.6f" % distribution[testCase].computeEntropy())

    # Confidence regions
    if distribution[testCase].getDimension() == 1:
        interval, threshold = distribution[
            testCase
        ].computeMinimumVolumeIntervalWithMarginalProbability(0.95)
        prec = ot.PlatformInfo.GetNumericalPrecision()
        ot.PlatformInfo.SetNumericalPrecision(5)
        print("Minimum volume interval=", interval)
        ot.PlatformInfo.SetNumericalPrecision(prec)
        print("threshold=", ot.Point(1, threshold))
        levelSet, beta = distribution[
            testCase
        ].computeMinimumVolumeLevelSetWithThreshold(0.95)
        print("Minimum volume level set=", levelSet)
        print("beta=", ot.Point(1, beta))
        interval, beta = distribution[
            testCase
        ].computeBilateralConfidenceIntervalWithMarginalProbability(0.95)
        print("Bilateral confidence interval=", interval)
        print("beta=", ot.Point(1, beta))
        interval, beta = distribution[
            testCase
        ].computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, False)
        print("Unilateral confidence interval (lower tail)=", interval)
        print("beta=", ot.Point(1, beta))
        interval, beta = distribution[
            testCase
        ].computeUnilateralConfidenceIntervalWithMarginalProbability(0.95, True)
        print("Unilateral confidence interval (upper tail)=", interval)
        print("beta=", ot.Point(1, beta))

    mean = distribution[testCase].getMean()
    print("mean      =", repr(mean))
    print("mean (ref)=", repr(referenceDistribution[testCase].getMean()))
    standardDeviation = distribution[testCase].getStandardDeviation()
    print("standard deviation      =", repr(standardDeviation))
    print(
        "standard deviation (ref)=",
        repr(referenceDistribution[testCase].getStandardDeviation()),
    )
    skewness = distribution[testCase].getSkewness()
    print("skewness      =", repr(skewness))
    print("skewness (ref)=", repr(referenceDistribution[testCase].getSkewness()))
    kurtosis = distribution[testCase].getKurtosis()
    print("kurtosis      =", repr(kurtosis))
    print("kurtosis (ref)=", repr(referenceDistribution[testCase].getKurtosis()))
    covariance = distribution[testCase].getCovariance()
    print("covariance      =", repr(cleanCovariance(covariance)))
    print(
        "covariance (ref)=",
        repr(cleanCovariance(referenceDistribution[testCase].getCovariance())),
    )
    parameters = distribution[testCase].getParametersCollection()
    print("parameters      =", repr(parameters))
    print(
        "parameters (ref)=",
        repr(referenceDistribution[testCase].getParametersCollection()),
    )
    print("parameter       =", repr(distribution[testCase].getParameter()))
    print("parameter desc  =", repr(distribution[testCase].getParameterDescription()))
    print("marginal 0      =", repr(distribution[testCase].getMarginal(0)))
    print(
        "Standard representative=",
        referenceDistribution[testCase].getStandardRepresentative(),
    )

    mean = distribution[testCase].getMean()
    cpdf = distribution[testCase].computeConditionalPDF(mean[0], [])
    print(f"conditional PDF={cpdf:.6f}")
    scpdf = distribution[testCase].computeSequentialConditionalPDF(mean)
    print(f"sequential conditional PDF={scpdf}")
    ccdf = distribution[testCase].computeConditionalCDF(mean[0], [])
    print(f"conditional CDF={ccdf:.6f}")
    sccdf = distribution[testCase].computeSequentialConditionalCDF(mean)
    print(f"sequential conditional CDF={sccdf}")

    ot.Log.Show(ot.Log.TRACE)
    ot.RandomGenerator.SetSeed(1)
    validation = ott.DistributionValidation(distribution[testCase])
    validation.skipMinimumVolumeLevelSet()  # slow
    validation.run()

# Check simplification
candidates = [
    ot.Normal(1.0, 2.0),
    ot.Uniform(1.0, 2.0),
    ot.Exponential(1.0, 2.0),
    ot.TruncatedDistribution(ot.WeibullMin(), 1.5, 7.8),
    ot.Beta(1.5, 6.3, -1.0, 2.0),
    ot.JointDistribution([ot.Normal()] * 2),
    ot.BlockIndependentDistribution([ot.Normal(2), ot.Normal(2)]),
    ot.BlockIndependentCopula([ot.NormalCopula(2), ot.NormalCopula(2)]),
    ot.Dirichlet([0.7, 0.3]),
]
intervals = [
    ot.Interval(-1.0, 4.0),
    ot.Interval(0.2, 2.4),
    ot.Interval(2.5, 65.0),
    ot.Interval(2.5, 6.0),
    ot.Interval(-2.5, 6.0),
    ot.Interval(2),
    ot.Interval(4),
    ot.Interval(4),
    ot.Interval(0.2, 2.4),
]
for i in range(len(candidates)):
    d = ot.TruncatedDistribution(candidates[i], intervals[i])
    print("d=", d, "simplified=", d.getSimplifiedVersion())

# Check that we can set the bounds independently
truncated = ot.TruncatedDistribution()
truncated.setDistribution(ot.Normal(20.0, 7.5))
truncated.setBounds(ot.Interval([0], [80], [True], [False]))
print("after setbounds q@0.9=", truncated.computeQuantile(0.9))
# Test for issue #1190
dist = ot.Normal(6.3e-19, 2.1e-19)
dist = ot.TruncatedDistribution(dist, 4.2e-19, ot.TruncatedDistribution.LOWER)

# non-finite bound bug
bounds = ot.Interval([-2, -3], [2, 3.0], [True, False], [True, True])
dist = ot.TruncatedDistribution(ot.Normal(2), bounds)
print("proba=%.6f" % dist.computeCDF([3.0, -3.0]))

# relative range wrt quantile epsilon issue
unif = ot.Uniform(0.0, 1e12)
trunc = ot.TruncatedDistribution(unif, 0.25, 2.0)
print("q@0.1=", trunc.computeQuantile(0.1))

# n-d CDF inversion
dim = 10
unif = ot.JointDistribution([ot.Uniform(0.0, 1.0)] * dim)
trunc = ot.TruncatedDistribution(unif, ot.Interval([0.0] * dim, [1e-6] * dim))
x = trunc.getRealization()

# marginal of truncated is not truncated marginal
R = ot.CorrelationMatrix(2)
R[1, 0] = 0.8
normal = ot.Normal([0.0] * 2, [1.0] * 2, R)
trunc = ot.TruncatedDistribution(normal, ot.Interval([-0.5] * 2, [1.0] * 2))
marginal0 = trunc.getMarginal(0)
ott.assert_almost_equal(marginal0.getMean(), [0.220527])
