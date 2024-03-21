.. _chaos_basis:

Polynomial chaos basis
----------------------

This page is focused on a specific kind of functional chaos
expansion, namely
*polynomial chaos expansion*.

Introduction
~~~~~~~~~~~~

We consider the notations introduced in :ref:`functional_chaos`. We recall here some of them to facilitate the reading of the documentation.
Let :math:`T: \Rset^{n_X} \rightarrow \Rset^{n_X}` be an iso-probabilistic such that:

.. math::

    \vect{Z} = T(\vect{X}) \sim \mu_{\vect{Z}}

where the measure :math:`\mu_{\vect{Z}}` is uniquely defined by all its moments and has independent components.
We consider the model :math:`h = g \circ T` and we assume that
:math:`h(\vect{Z}) \in \Rset`. However, the
following derivations hold in case of a multivariate output.
Let us assume that:

-  :math:`Y = h(\vect{Z})` has a finite variance, i.e.
   :math:`\Var{h(\vect{Z})} < \infty`;

-  :math:`\vect{Z}` has independent components.


Polynomial chaos expansion
~~~~~~~~~~~~~~~~~~~~~~~~~~

The meta model called *polynomial chaos expansion* is a particular *functional chaos expansion*
which makes the following choices:

- the projection spaces :math:`\cP_n` are a sequence of nested polynomial subspaces:
  :math:`\cP_n \subset \cP_{n+1}`. Now, for these polynomial spaces :math:`\cP_n` to be subspaces
  of :math:`L^2(\mu_{\vect{Z}})`, it is necessary that the measure :math:`\mu_{\vect{Z}}` has all its
  finite moments;

- the particular basis of :math:`\cP_n` consists of the family of orthonormal polynomials with respect
  to the measure :math:`\mu_{\vect{Z}}`. The basis of :math:`\cP_{n+1}` is constructed from that of
  :math:`\cP_n` by completing it.


To ensure :eq:`fermeturePn`, the measure :math:`\mu_{\vect{Z}}` needs to be uniquely defined by all its moments.

As :math:`\vect{Z}` has independent components, we have:

  .. math::

     \mu_{\vect{Z}}(\vect{z})= \prod_{i=1}^d \mu_i(z_i)

      
Thus, the multivariate polynomial basis can be built using the *tensor product* of the univariate polynomial bases which are orthonormal with respect to
the marginals of :math:`\mu_{\vect{Z}}`. Introducing the multi-index :math:`\vect{\alpha} = (\alpha_1, \dots, \alpha_d)` representing the marginal polynomial degrees, we built the basis of orthonormal polynomials with respect to :math:`\mu_{\vect{Z}}`
as follows:

  .. math::

        \Psi_\vect{\alpha}(\vect{z}) = \prod_{i=1}^d \Psi_{\alpha_i}(z_i)


The orthonormal polynomial basis with respect to the marginal :math:`\mu_i` is known for some distributions.
This is the case for:

- the measure :math:`\mu_i = \cU` whose orthonormal polynomial basis is the family of Legendre polynomials,

- the measure :math:`\mu_i = \Gamma` whose orthonormal polynomial basis is the family of Laguerre polynomials,

- the measure :math:`\mu_i = \beta` whose orthonormal polynomial basis is the family of Jacobi polynomials,

- the measure :math:`\mu_i = \cN` whose orthonormal polynomial basis is the family of Hermitte polynomials.

If the family is not already known, the polynomials can be represented by their three-term
recurrence, using the adaptive Stieljes algorithm (see :class:`~openturns.AdaptiveStieltjesAlgorithm`).
Once the sequence of recurrence coefficients
is known, the Reverse Clenshaw algorithm enables fast, stable evaluation of the polynomials
at any point (see :class:`~openturns.OrthogonalUniVariatePolynomial`).

Then, the meta model of *h* is the solution of:

  .. math::
    :label: PC

     \widetilde{h}   = \argmin_{h_n \in \cP_n} \| f-h_n \|^2_{L^2(\mu_{\vect{Z}})} = \sum_{k \in I_n}  a_k \Psi_k


where :math:`\{a_k \in \Rset\}_{k\in I_n}` are real coefficients.

Several strategies are possible to compute the coefficients :math:`(a_k)_k` (see :ref:`response_surface` to get more details on the computation of the coefficents):

- The orthonormality property of the polynomial basis makes it easy to calculate the coefficients which
are defined by scalar products (see :class:`~openturns.IntegrationExpansion`):

  .. math::
      :label: coeffAlphak

      a_k = \langle g,  \psi_k \rangle = \Expect{g(\vect{Z}) \psi_k(\vect{Z})}


- The coefficients are solution of the discretized least square problem :eq:`PC` (see :class:`~openturns.LeastSquaresExpansion`).

.. topic:: API:

    - See :class:`~openturns.AdaptiveStieltjesAlgorithm`
    - See :class:`~openturns.OrthogonalUniVariatePolynomial`
    - See :class:`~openturns.OrthogonalUniVariatePolynomialFactory`
    - See :class:`~openturns.OrthogonalUniVariatePolynomialFamily`
    - See :class:`~openturns.IntegrationExpansion`
    - See :class:`~openturns.LeastSquaresExpansion`


.. topic:: Examples:

    - See :doc:`/auto_meta_modeling/polynomial_chaos_metamodel/plot_functional_chaos`


.. topic:: References:

    - [soizeghanem2004]_
    - [ghanem1991]_
    - [lemaitre2010]_
    - R. Ghanem and P. Spanos, 1991, "Stochastic finite elements -- A spectral approach", Springer Verlag. (Reedited by Dover Publications, 2003).
