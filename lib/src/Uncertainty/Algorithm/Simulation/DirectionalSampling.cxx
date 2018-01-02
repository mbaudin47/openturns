//                                               -*- C++ -*-
/**
 *  @brief DirectionalSampling is an implementation of the directional
 *
 *  Copyright 2005-2018 Airbus-EDF-IMACS-Phimeca
 *
 *  This library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/DirectionalSampling.hxx"
#include "openturns/CorrelationMatrix.hxx"
#include "openturns/LinearFunction.hxx"
#include "openturns/ComposedFunction.hxx"
#include "openturns/Matrix.hxx"
#include "openturns/ComparisonOperatorImplementation.hxx"
#include "openturns/GaussLegendre.hxx"

BEGIN_NAMESPACE_OPENTURNS

/*
 * @class DirectionalSampling
 */

CLASSNAMEINIT(DirectionalSampling)

static const Factory<DirectionalSampling> Factory_DirectionalSampling;

/* Constructor with parameters */
DirectionalSampling::DirectionalSampling()
  : Simulation()
  , standardFunction_(standardEvent_.getImplementation()->getFunction())
  , inputDistribution_(standardEvent_.getImplementation()->getAntecedent()->getDistribution().getImplementation())
{
  // Nothing to do
}

/* Constructor with parameters */
DirectionalSampling::DirectionalSampling(const Event & event)
  : Simulation(event)
{
  if (!event.isComposite()) throw InvalidArgumentException(HERE) << "DirectionalSampling requires a composite event";
  standardEvent_ = StandardEvent(event);
  standardFunction_ = standardEvent_.getImplementation()->getFunction();
  inputDistribution_ = standardEvent_.getImplementation()->getAntecedent()->getDistribution().getImplementation();
  samplingStrategy_ = SamplingStrategy(inputDistribution_->getDimension());
}

/* Constructor with parameters */
DirectionalSampling::DirectionalSampling(const Event & event,
    const RootStrategy & rootStrategy,
    const SamplingStrategy & samplingStrategy)
  : Simulation(event)
  , rootStrategy_(rootStrategy)
{
  if (!event.isComposite()) throw InvalidArgumentException(HERE) << "DirectionalSampling requires a composite event";
  standardEvent_ = StandardEvent(event);
  standardFunction_ = standardEvent_.getImplementation()->getFunction();
  inputDistribution_ = standardEvent_.getImplementation()->getAntecedent()->getDistribution().getImplementation();
  setSamplingStrategy(samplingStrategy);
}

/* Virtual constructor */
DirectionalSampling * DirectionalSampling::clone() const
{
  return new DirectionalSampling(*this);
}

/* Compute the contribution of a direction to the probability given the roots x_0,...,x_{n-1} of the performance function along the direction.
   If the origin is in the failure space:
   dP = 1.0 - \sum_{k=0}^{n-1}F^c(x_k)
   If the origin is not in the failure space:
   dP = \sum_{k=0}^{n-1}F^c(x_k)
*/
Scalar DirectionalSampling::computeContribution(const ScalarCollection & roots)
{
  Scalar sign = 1.0;
  Scalar estimate = 0.0;
  const UnsignedInteger size = roots.getSize();
  for (UnsignedInteger indexRoot = 0; indexRoot < size; ++indexRoot)
  {
    Scalar currentRoot = roots[indexRoot];
    estimate += sign * inputDistribution_->computeRadialDistributionCDF(currentRoot, true);
    sign = -sign;
  }
  // Is the origin in the failure space?
  // Here, we know that the getOriginValue() method will not throw an exception, as we already called the solve() method
  // of the root strategy, which in turn initialized the computation of the origin value.
  if (standardEvent_.getDomain().contains(Point(1, rootStrategy_.getOriginValue()))) estimate = 1.0 - estimate;
  return estimate;
}

/* Compute the mean point of a direction given the roots x_0,...,x_{n-1} of the performance function along the direction.
   If the origin is in the failure space we add a root at 0, and if the resulting number of roots is odd we add a root at +\infty.
   The integrals \int_{x_k}^{x_{k+1}} xp(x)dx = -[xF^c(x)]_{x_k}^{x_{k+1}} + \int_{x_k}^{x_{k+1}} F^c(x)dx are computed using a Gauss-Legendre quadrature rule.
*/
Scalar DirectionalSampling::computeMeanContribution(const ScalarCollection & roots)
{
  ScalarCollection xK(0);
  // Is the origin in the failure space?
  // Here, we know that the getOriginValue() method will not throw an exception, as we already called the solve() method
  // of the root strategy, which in turn initialized the computation of the origin value.
  if (standardEvent_.getDomain().contains(Point(1, rootStrategy_.getOriginValue()))) xK.add(0.0);
  const UnsignedInteger size = roots.getSize();
  for (UnsignedInteger indexRoot = 0; indexRoot < size; ++indexRoot) xK.add(roots[indexRoot]);
  // If the number of points is odd, add a point at infinity
  if (xK.getSize() % 2 == 1) xK.add(rootStrategy_.getMaximumDistance());
  // Here we know that the number of points is even. We can integrate the contributions.
  const UnsignedInteger segmentNumber = xK.getSize() / 2;
  // Quadrature rule
  const UnsignedInteger integrationNodesNumber = ResourceMap::GetAsUnsignedInteger("DirectionalSampling-MeanContributionIntegrationNodesNumber");
  const GaussLegendre integrator(Indices(1, integrationNodesNumber));
  const Point nodes(integrator.getNodes().getImplementation()->getData() * 2.0 - Point(integrationNodesNumber, 1.0));
  const Point weights(integrator.getWeights() * 2.0);
  Scalar value = 0.0;
  for (UnsignedInteger segmentIndex = 0; segmentIndex < segmentNumber; ++segmentIndex)
  {
    const Scalar a = xK[2 * segmentIndex];
    const Scalar b = xK[2 * segmentIndex + 1];
    const Scalar halfLength = 0.5 * (b - a);
    // Accumulate the bracket part
    value += a * inputDistribution_->computeRadialDistributionCDF(a, true) - b * inputDistribution_->computeRadialDistributionCDF(b, true);
    // Compute the integral part
    Scalar sum = 0.0;
    for (UnsignedInteger k = 0; k < integrationNodesNumber; ++k) sum += weights[k] * inputDistribution_->computeRadialDistributionCDF(a + (1.0 + nodes[k]) * halfLength, true);
    sum *= halfLength;
    // Accumulate the integral part
    value += sum;
  }
  return value;
}

/* Compute the contribution of a set of directions direction to the probability */
Scalar DirectionalSampling::computeTotalContribution(const Sample & directionSample)
{
  const UnsignedInteger sampleSize = directionSample.getSize();
  const UnsignedInteger dimension = directionSample.getDimension();
  Scalar totalContribution = 0.0;
  // meanPointInEventDomain = Point(dimension);
  UnsignedInteger contributionNumber = 0;
  Matrix linear(dimension, 1);
  // For each direction
  for (UnsignedInteger indexDirection = 0; indexDirection < sampleSize; ++indexDirection)
  {
    const Point direction(directionSample[indexDirection]);
    // First compute the roots along this direction
    // 1. Build the scalar function along the direction
    // 1.1 Build the linear function along the direction
    for (UnsignedInteger indexComponent = 0; indexComponent < dimension; ++indexComponent) linear(indexComponent, 0) = direction[indexComponent];
    const LinearFunction ray(Point(1, 0.0), Point(dimension, 0.0), linear);
    // 1.2 Build the function along the ray
    const ComposedFunction functionAlongRay(standardFunction_, ray);
    // 2. Solve the function along the ray
    const ScalarCollection roots(rootStrategy_.solve(functionAlongRay, standardEvent_.getThreshold()));
    // Second, compute the contribution of this direction
    const Scalar contribution = computeContribution(roots);
    // If there is a contribution in this direction
    if (contribution > 0.0)
    {
      // Accumulate the contribution
      totalContribution += contribution;
      // Third, compute the mean point along this direction
      // meanPointInEventDomain = meanPointInEventDomain + computeMeanContribution(roots) * direction;
      ++contributionNumber;
    } // if contribution
  }
  // meanPointInEventDomain = meanPointInEventDomain * (1.0 / contributionNumber);
  return totalContribution / sampleSize;
}

/* Compute the block sample and the points that realized the event */
Sample DirectionalSampling::computeBlockSample()
{
  const UnsignedInteger size = getBlockSize();
  Sample blockSample(size, 1);
  // For each entry of the block sample
  // realizedEventSample = Sample(blockSize_, event_.getImplementation()->getAntecedent()->getDistribution().getDimension());
  for (UnsignedInteger index = 0; index < size; ++index)
  {
    const Sample directionSample(samplingStrategy_.generate());
    // Compute the contribution of the sub-sample computed according to the sampling strategy
    const Scalar contribution = computeTotalContribution(directionSample);
    blockSample[index][0] = contribution;
  }
  return blockSample;
}

/* Root strategy accessor */
void DirectionalSampling::setRootStrategy(const RootStrategy & rootStrategy)
{
  rootStrategy_ = rootStrategy;
}

RootStrategy DirectionalSampling::getRootStrategy() const
{
  return rootStrategy_;
}

/* Sampling strategy */
void DirectionalSampling::setSamplingStrategy(const SamplingStrategy & samplingStrategy)
{
  samplingStrategy_ = samplingStrategy;

  // To force the sampling strategy to have the correct dimension
  samplingStrategy_.setDimension(inputDistribution_->getDimension());
}

SamplingStrategy DirectionalSampling::getSamplingStrategy() const
{
  return samplingStrategy_;
}

/* String converter */
String DirectionalSampling::__repr__() const
{
  OSS oss;
  oss << "class=" << DirectionalSampling::GetClassName()
      << " rootStrategy=" << rootStrategy_.__repr__()
      << " samplingStrategy=" << samplingStrategy_.__repr__()
      << " derived from " << Simulation::__repr__();
  return oss;
}

/* Method save() stores the object through the StorageManager */
void DirectionalSampling::save(Advocate & adv) const
{
  Simulation::save(adv);
  adv.saveAttribute("rootStrategy_", rootStrategy_);
  adv.saveAttribute("samplingStrategy_", samplingStrategy_);
}

/* Method load() reloads the object from the StorageManager */
void DirectionalSampling::load(Advocate & adv)
{
  Simulation::load(adv);
  adv.loadAttribute("rootStrategy_", rootStrategy_);
  adv.loadAttribute("samplingStrategy_", samplingStrategy_);
  standardEvent_ = StandardEvent(event_);
  standardFunction_ = standardEvent_.getImplementation()->getFunction();
  inputDistribution_ = standardEvent_.getImplementation()->getAntecedent()->getDistribution().getImplementation();
}

END_NAMESPACE_OPENTURNS
