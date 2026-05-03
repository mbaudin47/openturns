//                                               -*- C++ -*-
/**
 *  @brief Quantile confidence interval computation
 *
 *  Copyright 2005-2026 Airbus-EDF-IMACS-ONERA-Phimeca
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
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <iomanip>
#include "openturns/QuantileConfidence.hxx"
#include "openturns/PersistentObjectFactory.hxx"
#include "openturns/Binomial.hxx"
#include "openturns/SpecFunc.hxx"
#include "openturns/Brent.hxx"
#include "openturns/SymbolicFunction.hxx"
#include "openturns/DistFunc.hxx"
#include "openturns/ParametricFunction.hxx"
#include "openturns/EvaluationImplementation.hxx"

BEGIN_NAMESPACE_OPENTURNS



CLASSNAMEINIT(QuantileConfidence)

static const Factory<QuantileConfidence> Factory_QuantileConfidence;

/*
 * @class QuantileConfidence
*/

QuantileConfidence::QuantileConfidence()
  : PersistentObject()
{
  // Nothing to to do
}

/* Constructor with parameters */
QuantileConfidence::QuantileConfidence(const Scalar alpha, const Scalar beta)
  : PersistentObject()
  , alpha_(alpha)
  , beta_(beta)
{
  if (!(alpha >= 0.0) || !(alpha <= 1.0))
    throw InvalidArgumentException(HERE) << "Quantile level must be in [0, 1]";
  if (!(beta >= 0.0) || !(beta <= 1.0))
    throw InvalidArgumentException(HERE) << "Confidence level must be in [0, 1]";
}

/* Virtual constructor */
QuantileConfidence * QuantileConfidence::clone() const
{
  return new QuantileConfidence(*this);
}

/* String converter */
String QuantileConfidence::__repr__() const
{
  return OSS() << GetClassName()
         << " alpha=" << alpha_
         << " beta=" << beta_;
}

String QuantileConfidence::__str__(const String & /*offset*/) const
{
  return OSS() << GetClassName() << "(alpha=" << alpha_ << ", beta=" << beta_ << ")";
}

// compute rank k given by the beta-level quantile of B(n, alpha)
UnsignedInteger QuantileConfidence::computeUnilateralRank(const UnsignedInteger size, const Bool lower_bounded) const
{
  const UnsignedInteger minimumSize = computeUnilateralMinimumSampleSize(0, lower_bounded);
  if (size < minimumSize)
    throw InvalidArgumentException(HERE) << "Cannot compute unilateral rank as size (" << size << ") is lower than minimum size (" << minimumSize << ")";

  const Binomial binomial(size, alpha_);
  const Scalar p = lower_bounded ? 1.0 - std::pow(1.0 - alpha_, 1.0 * size) : 1.0 - std::pow(alpha_, 1.0 * size);

  if (p < beta_)
  {
    if (lower_bounded)
      throw InvalidArgumentException(HERE) << "Cannot compute rank as parameters do not satisfy 1 - alpha^n >= beta";
    else
      throw InvalidArgumentException(HERE) << "Cannot compute rank as parameters do not satisfy 1 - (1 - alpha)^n >= beta";
  }
  const UnsignedInteger rank = binomial.computeQuantile(beta_, lower_bounded)[0];
  return rank;
}

// compute argmin (k1, k2) P_X([k1, k2]) under constraint P_X([k1, k2])>=beta, with X~B(n, alpha)
Indices QuantileConfidence::computeBilateralRank(const UnsignedInteger size) const
{
  const UnsignedInteger minimumSize = computeBilateralMinimumSampleSize();
  if (size < minimumSize)
    throw InvalidArgumentException(HERE) << "Cannot compute bilateral rank as size (" << size << ") is lower than minimum size (" << minimumSize << ")";

  // Consider the Binomial(n=size, alpha).
  // Find the indices (k1, k2) with smallest probability >=beta.
  const Binomial binomial(size, alpha_);
  Scalar pBest = 1.0;
  UnsignedInteger k1Best = 0;
  UnsignedInteger k2Best = size;
  UnsignedInteger k1 = 0;
  Bool found = false;
  while (k1 < size)
  {
    // Compute p1 = CDF(k)
    const Scalar p1 = binomial.computeCDF(k1);

    // Check that CDF(k) + beta <= 1, otherwise computeScalarQuantile fail
    if (p1 + beta_ > 1.0)
      break;

    // Compute k2 from p1
    const UnsignedInteger k2 = binomial.computeScalarQuantile(p1 + beta_);

    // Compute P(k1 < X <= k2) = P(k1 + 1 <= X <= k2)
    // Here, CDF(k2) - CDF(k1) also works
    const Scalar p = binomial.computeProbability(Interval(k1 + 1, k2));
    if ((p >= beta_) && (p < pBest))
    {
      pBest = p;
      k1Best = k1;
      k2Best = k2;
      found = true;
    }
    // Find the next significant probability jump
    UnsignedInteger ell;
    const Bool succeeds = searchProbabilityJump(size, k1, ell);
    if (succeeds)
      k1 = ell;
    else
      break;
    
  } // while k1
  if (!found)
    throw InvalidArgumentException(HERE) << "Cannot find suitable ranks for size=" 
        << size << ", alpha=" << alpha_ << ", beta=" << beta_;
  return {k1Best, k2Best};
}

/** Find the smallest ell such that F(ell) > F(k) + epsilon and ell > k */
Bool QuantileConfidence::searchProbabilityJump(const UnsignedInteger size, const UnsignedInteger k, UnsignedInteger & ell) const
{
  const Scalar epsilon_p = ResourceMap::GetAsScalar("QuantileConfidence-ProbabilityEpsilon");
  const Binomial binomial(size, alpha_);
  const Scalar p_ref = binomial.computeCDF(k) + epsilon_p;
  Bool succeeds = false;

  if (p_ref >= 1.0)
  {
    ell = size;
    return succeeds;
  }

  // Phase 1: Increment increase following the 2^i law
  UnsignedInteger i = 0;
  UnsignedInteger increment = 1;

  while (true)
  {
    const UnsignedInteger L_test = std::min(k + increment, size);
    const Scalar p_test = binomial.computeCDF(L_test);

    if (p_test > p_ref)
      break;

    if (L_test == size)
    {
      ell = size;
      return succeeds;
    }

    ++i;
    increment *= 2;
  }

  // Phase 2: Increment decrease following a law of the form 2^-j
  UnsignedInteger L_base = (i > 0) ? k + (increment / 2) : k;
  UnsignedInteger reduction = increment;

  while (true)
  {
    reduction /= 2;
    if (reduction == 0)
      break;

    const UnsignedInteger L_milieu = L_base + reduction;
    const Scalar p_milieu = binomial.computeCDF(L_milieu);

    if (p_milieu <= p_ref)
    {
      // If no jump is detected, the starting point is moved forward
      L_base = L_milieu;
    }
  }

  succeeds = true;
  ell = L_base + 1;
  return succeeds;
}

// compute interval of the form [X_k; +inf[ or ]-inf; X_k] from unilateral rank k
Interval QuantileConfidence::computeUnilateralConfidenceInterval(const Sample & sample, const Bool lower_bounded) const
{
  Scalar coverageOut = -1.0;
  return computeUnilateralConfidenceIntervalWithCoverage(sample, coverageOut, lower_bounded);
}

// compute interval of the form [X_k; +inf[ or ]-inf; X_k] from unilateral rank k and the actual coverage
Interval QuantileConfidence::computeUnilateralConfidenceIntervalWithCoverage(const Sample & sample, Scalar & coverageOut, const Bool lower_bounded) const
{
  if (sample.getDimension() != 1)
    throw InvalidArgumentException(HERE) << "Expected a sample of dimension 1";
  Point lowerBound(1, -SpecFunc::MaxScalar);
  Point upperBound(1, SpecFunc::MaxScalar);
  Interval::BoolCollection finiteLowerBound(1, false);
  Interval::BoolCollection finiteUpperBound(1, false);
  const Sample sortedSample(sample.sort());
  const UnsignedInteger k = computeUnilateralRank(sample.getSize(), lower_bounded);
  if (lower_bounded)
  {
    lowerBound[0] = sortedSample(k, 0);
    finiteLowerBound[0] = true;
  }
  else
  {
    upperBound[0] = sortedSample(k, 0);
    finiteUpperBound[0] = true;
  }
  coverageOut = computeUnilateralCoverage(sample.getSize(), k, lower_bounded);
  return Interval(lowerBound, upperBound, finiteLowerBound, finiteUpperBound);
}

// compute interval of the form [X_k1; X_k2] from bilateral ranks k1, k2
Interval QuantileConfidence::computeBilateralConfidenceInterval(const Sample & sample) const
{
  Scalar coverageOut = -1.0;
  return computeBilateralConfidenceIntervalWithCoverage(sample, coverageOut);
}

// compute interval of the form [X_k1; X_k2] from bilateral ranks k1, k2 with actual coverage
Interval QuantileConfidence::computeBilateralConfidenceIntervalWithCoverage(const Sample & sample, Scalar & coverageOut) const
{
  if (sample.getDimension() != 1)
    throw InvalidArgumentException(HERE) << "Expected a sample of dimension 1";
  const Indices rank(computeBilateralRank(sample.getSize()));
  coverageOut = computeBilateralCoverage(sample.getSize(), rank[0], rank[1]);
  const Sample sortedSample(sample.sort());
  return Interval(Point({sortedSample(rank[0], 0)}), Point({sortedSample(rank[1], 0)}));
}

Scalar QuantileConfidence::computeUnilateralCoverage(const UnsignedInteger size, const UnsignedInteger rank, const Bool lower_bounded) const
{
  if (rank >= size)
    throw InvalidArgumentException(HERE) << "The rank must be strictly less than the sample size.";

  // Y ~ Binomial(size, alpha)
  Binomial binomial(size, alpha_);  
  Scalar cdf;

  if (lower_bounded)
    // Lower bound coverage: 1 - P(Y <= rank)
    cdf = binomial.computeComplementaryCDF(rank);
  else
    // Upper bound coverage: P(Y <= rank)
    cdf = binomial.computeCDF(rank);
  
  return cdf;
}

Scalar QuantileConfidence::computeBilateralCoverage(const UnsignedInteger size, const UnsignedInteger rank1, const UnsignedInteger rank2) const
{
  if (rank1 > rank2)
    throw InvalidArgumentException(HERE) << "The lower rank (" << rank1 
                                         << ") must be strictly less than the upper rank (" << rank2 << ").";
  if (rank2 >= size)
    throw InvalidArgumentException(HERE) << "The upper rank (" << rank2 
                                         << ") must be strictly less than the sample size (" << size << ").";

  // Y ~ Binomial(size, alpha)
  Binomial binomial(size, alpha_);
  
  // P(k1 + 1 <= Y <= k2) = P(Y <= k2) - P(Y <= k1)
  Scalar cdf_upper = binomial.computeCDF(rank2);
  Scalar cdf_lower = binomial.computeCDF(rank1);
  
  return cdf_upper - cdf_lower;
}

namespace
{
class QuantileConfidenceEvaluation: public EvaluationImplementation
{
public:
  QuantileConfidenceEvaluation(const Scalar alpha,
                               const UnsignedInteger rank,
                               const Bool lower_bounded)
    : EvaluationImplementation()
    , alpha_(alpha)
    , rank_(rank)
    , tail_(lower_bounded)
  {
    // Nothing to do
  }

  QuantileConfidenceEvaluation * clone() const override
  {
    return new QuantileConfidenceEvaluation(*this);
  }

  Point operator() (const Point & point) const override
  {
    return {DistFunc::pBeta(rank_ + 1, point[0] - rank_, tail_ ? alpha_ : 1.0 - alpha_)};
  }

  UnsignedInteger getInputDimension() const override
  {
    return 1;
  }

  UnsignedInteger getOutputDimension() const override
  {
    return 1;
  }

private:
  Scalar alpha_ = 0.0;
  UnsignedInteger rank_ = 0;
  Bool tail_ = true;
}; // class QuantileConfidenceEvaluation
} // Anonymous namespace


UnsignedInteger QuantileConfidence::computeUnilateralMinimumSampleSize(const UnsignedInteger tail_rank, const Bool lower_bounded) const
{
  // Here we have to find the minimal value of N such that
  // 1-\sum_{i=N-r}^N Binomial(i, N)\alpha^i(1-\alpha)^{N-i}>=\beta
  // where:
  //   \alpha = alpha_
  //   \beta  = beta_
  //   r      = rank
  // It rewrites F_{N,alpha}(N-r-1)>=beta where F_{N,alpha} is the CDF of the
  // Binomial(N, alpha) distribution.
  // Easy case: rank=0, the quantile bound is given by the largest upper statistics. The equation to solve is N=\min{n|1-\alpha^n>=\beta}
  Scalar nApprox = 0.0;
  const Function wilksConstraint(QuantileConfidenceEvaluation(alpha_, tail_rank, lower_bounded).clone());
  if (tail_rank == 0)
  {
    if (lower_bounded)
      nApprox = std::log1p(-beta_) / std::log1p(-alpha_);
    else
      nApprox = std::log1p(-beta_) / std::log(alpha_);
  }
  else
  {
    // Search for the minimal sample size
    // First, search for an upper bound
    // Here we use the relation F_{N,alpha}(N - r - 1) = pBeta(N - k, k + 1, 1 - alpha)
    // where k = N - r - 1
    // We compute a reasonable guess for n using a Normal approximation:
    // n*alpha+q_beta*sqrt(n*alpha*(1-alpha))=n-r:
    const Scalar aBeta = DistFunc::qNormal(beta_);
    UnsignedInteger nMax = static_cast<UnsignedInteger>((tail_rank + 0.5 * (alpha_ * aBeta * aBeta + std::abs(aBeta) * std::sqrt(alpha_ * (4.0 * tail_rank + alpha_ * aBeta * aBeta)))) / (1.0 - alpha_));
    // This loop must end as wilksConstraint->1 when n->inf
    while (wilksConstraint(Point({0.0 + nMax}))[0] < beta_)
      // At the beginning of the loop nMax is >= 1 so it increases strictly
      nMax += nMax;
    nApprox = Brent().solve(wilksConstraint, beta_, tail_rank, nMax);
  } // rank > 0
  // Here, nApprox can be very close to an integer (in the sense of: the value of the constraint evaluated at the nearest integer is very close to beta_), in which case the correct answer is round(nApprox), or the answer is the next integer value
  const UnsignedInteger nInf = std::round(std::max(1.0 * tail_rank, nApprox));
  const Scalar constraintInf = wilksConstraint(Point({0.0 + nInf}))[0];
  if (std::abs(constraintInf - beta_) < std::sqrt(SpecFunc::Precision)) return nInf;
  return std::ceil(nApprox);
}

// Find the minimal value of N such that 1 - alpha^n - (1 - alpha)^n >= beta
UnsignedInteger QuantileConfidence::computeBilateralMinimumSampleSize() const
{
  const Scalar gamma = std::max(alpha_, 1.0 - alpha_);
  const Scalar nMin = std::ceil(std::log1p(-beta_) / std::log(gamma));
  const Scalar nMax = std::ceil((std::log1p(-beta_) - std::log(2)) / std::log(gamma));
  const SymbolicFunction residualFunction(Description({"n", "alpha", "beta"}), Description({"1 - alpha^n - (1 - alpha)^n - beta"}));
  const ParametricFunction residualParametric(residualFunction, Indices({1, 2}), Point({alpha_, beta_}));
  Brent solver;
  const Scalar root = solver.solve(residualParametric, 0.0, nMin - 1, nMax);
  return std::ceil(root);
}


// Compute k1, k2 such that lim n=>inf P(X_k1<x_alpha<X_k2)=beta, see [delmas2006] proposition 12.2.13 page 257
Indices QuantileConfidence::computeAsymptoticBilateralRank(const UnsignedInteger size) const
{
  const Scalar z = DistFunc::qNormal((1.0 + beta_) * 0.5);
  const Scalar delta = z * std::sqrt(alpha_ * (1.0 - alpha_) * size);
  const UnsignedInteger k1 = std::max(0.0, std::floor(size * alpha_ - delta - 1.0));
  const UnsignedInteger k2 = std::min(size - 1.0, std::floor(size * alpha_ + delta - 1.0));
  return {k1, k2};
}

// compute interval of the form [X_k1; X_k2] from asymptotic bilateral ranks k1, k2
Interval QuantileConfidence::computeAsymptoticBilateralConfidenceInterval(const Sample & sample) const
{
  const Indices rank(computeAsymptoticBilateralRank(sample.getSize()));
  const Sample sortedSample(sample.sort());
  return Interval(Point({sortedSample(rank[0], 0)}), Point({sortedSample(rank[1], 0)}));
}

/* Quantile level accessor */
void QuantileConfidence::setAlpha(const Scalar alpha)
{
  alpha_ = alpha;
}

Scalar QuantileConfidence::getAlpha() const
{
  return alpha_;
}

/* Confidence level accessor */
void QuantileConfidence::setBeta(const Scalar beta)
{
  beta_ = beta;
}

Scalar QuantileConfidence::getBeta() const
{
  return beta_;
}


/* Method save() stores the object through the StorageManager */
void QuantileConfidence::save(Advocate & adv) const
{
  PersistentObject::save(adv);
  adv.saveAttribute( "alpha_", alpha_ );
  adv.saveAttribute( "beta_", beta_ );
}

/* Method load() reloads the object from the StorageManager */
void QuantileConfidence::load(Advocate & adv)
{
  PersistentObject::load(adv);
  adv.loadAttribute( "alpha_", alpha_ );
  adv.loadAttribute( "beta_", beta_ );
}

END_NAMESPACE_OPENTURNS
