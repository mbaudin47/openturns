.. _functional_chaos:

Functional Chaos Expansion
--------------------------

Introduction
~~~~~~~~~~~~

Accounting for the joint probability density function (PDF)
:math:`\mu_{\vect{X}}(\vect{x})` of the input random vector
:math:`\vect{X}`, one seeks the joint PDF of output random vector
:math:`\vect{Y} = g(\vect{X})`. This may be achieved using
Monte Carlo (MC) simulation (see :ref:`monte_carlo_simulation`). However, the MC
method may require a large number of model evaluations, i.e. a great
computational cost, in order to obtain accurate results.

A possible solution to overcome this problem is to project the model
:math:`g` in a suitable functional space, such as
the Hilbert space :math:`L^2(\mu_{\vect{X}})` of square-integrable functions with
respect to the PDF :math:`\mu_{\vect{X}}`.
More precisely, we may consider an expansion of the model onto an orthonormal basis of :math:`L^2(\mu_{\vect{X}})`.
As an example of this type of expansion, one can mention expansions by
wavelets, polynomials, etc.

The principles of the building of a functional chaos expansion are described in the sequel.

Model
~~~~~

We consider the output random vector:

.. math::

    \vect{Y} = g(\vect{X})

where :math:`g: \Rset^{n_X} \rightarrow \Rset^{n_Y}` is the model,
:math:`\vect{X}` is the input random vector which distribution is
:math:`\mu_{\vect{X}}`,
:math:`n_X \in \Nset` is the input dimension,
:math:`n_Y \in \Nset` is the output dimension.
We assume that :math:`\vect{Y}` has finite variance i.e.
:math:`g\in :math:`L^2(\mu_{\vect{X}})`.

When :math:`n_Y > 1`, the functional chaos algorithm is used on each marginal
of :math:`\vect{Y}`, using the same multivariate orthonormal basis for
all the marginals.
Thus, the method is detailed here for a scalar output :math:`Y` and
:math:`g: \Rset^{n_X} \rightarrow \Rset`.

Iso-probabilistic transformation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`T: \Rset^{n_X} \rightarrow \Rset^{n_X}` be an isoprobabilistic transformation
such that :math:`\vect{Z} = T(\vect{X})` follows the measure :math:`\mu_{\vect{Z}}`.
Let :math:`h` be the function defined by the equation:

.. math::
    h = g \circ T^{-1}.

Therefore :math:`h \in L^2(\mu_{\vect{Z}})`.


Orthonormal basis with respect to a measure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The functional space :math:`L^2(\mu_{\vect{Z}})` is a Hilbert space  which admits a hilbertian basis
denoted by :math:`(\Psi_k)_{k \geq 0}`. Each function :math:`h \in L^2(\mu_{\vect{Z}})`
can be written as (see [lemaitre2010]_ page 39):

  .. math::
     :label: fctExph

      h = \sum_{k=0}^{+\infty} a_k \Psi_k

The measure :math:`\mu_{\vect{Z}}` defines in that space the following scalar product written for continuous and discrete variables:

  .. math::

        \forall (f,g) \in L^2(\mu_{\vect{Z}}), \quad \langle f, g \rangle _{L^2(\mu_{\vect{Z}})} = \Expect{f(\vect{X})g(\vect{X})} & =  \int f(\vect{x}) g(\vect{x})\, \mu_{\vect{Z}}(\vect{z}) d\vect{x} \\
        & = \sum_\vect{z} f(\vect{z}) g(\vect{z})\, \Prob{\vect{Z} = \vect{z}}

and the norm is defined by:

  .. math::

        \|f\|^2_{L^2(p_{\vect{X}})} = \Expect{\left[f(\vect{X})\right]^2} & = \int [f(\vect{z})]^2\, \mu_{\vect{Z}}(\vect{z}) d\vect{z} \\
            & = \sum_\vect{z} f^2(\vect{z}) \,\Prob{\vect{Z} = \vect{z}}

The basis :math:`(\Psi_k)_{k \in \Nset}` is orthonormal with respect to :math:`\mu_{\vect{Z}}` if it
verifies the following properties:

.. math::
   :label: orthonorm

    \langle \Psi_i, \Psi_{j}\rangle  = \int \Psi_i\vect{z} \Psi_{j}\vect{z}d\vect{z} = \delta_{i,j}


or:

.. math::
   :label: orthonormDisc

    \langle \Psi_i, \Psi_{j}\rangle  = \sum_\vect{z} = \delta_{i,j}


where :math:`\delta_{i,j} =1` is the Kronecker symbol:

.. math::

  \delta_{i,j}
  =
  \begin{cases}
  1 & \textrm{ if } i = j \\
  0 & \textrm{otherwise.}
  \end{cases}

See :ref:`orthogonal basis <orthogonal_basis>` to know the available orthonormal basis.

Functional chaos expansion
~~~~~~~~~~~~~~~~~~~~~~~~~~
Building a functional chaos expansion of :math:`h` consists in making the following choices:

- choice of a finite projection space :math:`\cP_n` such that:

  .. math::
       :label: fermeturePn

       \overline{\cup_{n\in \mathbb{N}} \cP_n} = L^2(\mu_{\vect{Z}})

  For example, we can choose the polynomials of total degree less than :math:`n`.

- choice of a basis of :math:`\cP_n` denoted by  :math:`(\Psi_k)_{k \in I_n}`:

  .. math::
       :label: Pn

       \cP_n = \mbox{span} (\Psi_k)_{k \in I_n}


  where :math:`I_n` is finite. Thus each element :math:`h_n \in\cP_n` can be written as:

  .. math::

     h_n = \sum_{k \in I_n} a_k \psi_k


  For example, we can choose the canonical basis or the family of orthonormal polynomials with respect to :math:`\mu_{\vect{Z}}`.


Then, the meta model of *h* is the solution of:

  .. math::
    :label: metaModeleh

     \widetilde{h}  = \argmin_{h_n \in \cP_n} \| f-h_n \|^2_{L^2(\mu_{\vect{Z}})}

which is a least-squares otimization problem.



The choice of the projection space :math:`\cP_n` and its basis :math:`(\Psi_k)_{k \in I_n}` is
designed to ensure that the discretized problem :eq:`metaModeleh` is easy to solve (well-conditioned
discrete problem).
In particular, the choice of basis has a major influence on the
conditioning of the least-squares problem :eq:`metaModeleh`.

Thus :math:`\widetilde{h}` is represented by a *finite* subset of coefficients :math:`(a_k)_{k\in I_n}` in a *truncated* basis :math:`(\Psi_k)_{k\in I_n}`:

.. math::

    \widetilde{h} = \sum_{k \in I_n}  a_k \Psi_k

The determination of :math:`I_n` can be made using one enumeration rule,
as presented in :ref:`enumeration_strategy`.
If the number of coefficients in :math:`I_n` is too large,
this can lead to *overfitting*.
This may happen e.g. if the total polynomial order we choose is too large.
In order to limit this effect, one method is to select the coefficients which
best predict the output, as presented in :ref:`polynomial_sparse_least_squares`.


**In OpenTURNS**, we choose a basis :math:`(\Psi_k)_{k \in I_n}` which is orthonormal with
respect to :math:`\mu_{\vect{Z}}`, so we have :eq:`orthonorm`. Furthermore, we require that the
first element be:

  .. math::
    :label: defPsi0

      \Psi_0 = 1

As for non-zero :math:`i`, :math:`\langle \psi_{i},\psi_{0} \rangle_{L^2(\mu_{\vect{Z}})} = 0`
by orthogonality of the base, relation :eq:`defPsi0` implies in particular that:

  .. math::

       \Expect{\psi_{i}(\vect{Z})} = \Expect{\Psi_{i}(\vect{Z})\Psi_{0}(\vect{Z})}= 0\quad \forall i\neq 0

The use of a basis orthonormal with respect to the measure :math:`\mu_{\vect{Z}}` facilitates the
computation of the :math:`a_k` coefficients, transforming the least-squares problem into a scalar product
calculation. In this case, the least-squares problem is equivalent to the computation of scalar products.
The algorithmic cost of solving the problem is much lower.


The meta model :math:`\widetilde{h}` can be used to build an efficient
random generator of :math:`Y` based on the random vector :math:`\vect{Z}`,
using the equation:

.. math::

    \widetilde{Y} = \widetilde{h}(\vect{Z})

This equation can be used to simulate independent random observations
from the functional chaos expansion.
This can be done by first simulating independent observations from
the distribution of the random vector :math:`\vect{Z}`,
then push forward these observations through the expansion.
See the :class:`~openturns.FunctionalChaosRandomVector` class
for more details on this topic.

Then, the meta model of *g* can be defined using the isoprobabilistic transformation :math:`T`:

.. math::
    :label: metaModeleg

    \widetilde{g} = \widetilde{h} \circ T

see  :ref:`response_surface` to get more details on:

- the available constructions of the truncated multivariate orthogonal basis,

- the computation of the coefficients.


Polynomial chaos expansion for independent variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OpenTURNS enables one to build the meta model called *polynomial chaos expansion* which makes the
following choices:

- the projection spaces :math:`\cP_n` are a sequence of nested polynomial subspaces:
  :math:`\cP_n \subset \cP_{n+1}`,

- the particular basis of :math:`\cP_n` consists of the family of orthonormal polynomials with respect
  to the measure :math:`\mu_{\vect{Z}}` if :math:`\mu_{\vect{Z}}` is such that the infinite sequel of
  its moments is defined.

Furthermore, to ensure :eq:`fermeturePn`, the measure :math:`\mu_{\vect{Z}}` needs to be uniquely defined
by all its moments. So, we poceed as follows:

- if the measure :math:`\mu_{\vect{X}}` is uniquely defined by all its moments, we use :math:`T=Id(\Rset^{n_X})`,

- if not, we use an iso-probabilistic transformation :math:`T` such that:

  .. math::
     :label: measureMu

     \vect{Z} = T(\vect{X})

is a random vector distributed according to the measure :math:`\mu_{\vect{Z}}` which is uniquely defined
by all its moments.
We also recommend to define :math:`\mu_{\vect{Z}}` with independent components in order to facilitate
the creation of the orthonormal basis as the tensorization of univariate polynomial basis orthonormal with
respect to its margins :math:`\mu_i` (see  :ref:`Polynomial chaos basis <chaos_basis>` and the classes
:class:`~openturns.OrthogonalUniVariatePolynomialFamily` and
:class:`~openturns.OrthogonalUniVariatePolynomialFactory`):

  .. math::

     \mu_{\vect{Z}}(\vect{z})= \prod_{i=1}^{n_X} \mu_i(z_i)


Other chaos expansions for independent variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After having mapped the input random vector :math:`\vect{X}` into the random vector :math:`\vect{Z}`
with independent components using  :math:`T` defined in :eq:`measureMu`, OpenTURNS enables one to use
the Haar wavelet functions or the Fourier series as orthonormal basis with respect to each margin
:math:`\mu_i`.

The Haar wavelets basis is orthonormal with respect to the the :math:`\cU(0,1)` measure (see
:class:`~openturns.HaarWaveletFactory`) and the Fourier series basis is orthonormal with respect to
the :math:`\cU(-\pi, \pi)` measure (see :class:`~openturns.FourierSeriesFactory`).


Some chaos expansions for dependent variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the components of the input random vector :math:`\vect{X}` are not independent, we can use an
iso-probabilistic transformation to map :math:`\vect{X}` into :math:`\vect{Z}` with independent components.

It is also possible to build up a multivariate orthonormal basis with respect to the
:math:`\mu_{\vect{X}}`  if it is uniquely defined by all its moments, as follows:

  .. math::

      \Psi_{\idx}(\vect{x}) \, \, = \,\,  K(\vect{x}) \;\prod_{i=1}^M \pi^{(i)}_{\alpha_{i}}(x_{i})


where :math:`K(\vect{x})` is a function of the copula of :math:`\vect{X}` and
:math:`\vect{\alpha} = (\alpha_1, \dots, \alpha_d)` a multi-index used to define the mutlivariate
polynomial basis built as the tensorization of the univariate orthonormal polynomial basis with
respect to :math:`\mu_i`  as follows:

  .. math::

        \Psi_\vect{\alpha}(\vect{z}) = \prod_{i=1}^d \Psi_{\alpha_i}(z_i).


OpenTURNS enables one to use the following kernel:

  .. math::

     K(\vect{x}) = \dfrac{1}{\sqrt{c(\vect{x}}}


where :math:`c` is the density of the copula of :math:`\vect{X}`. Then the orthonormal basis is
called the `Soize-Ghanem` basis (see
:class:`~openturns.SoizeGhanemFactory`).


Link with classical deterministic polynomial approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a deterministic setting (i.e. when the input parameters are
considered to be deterministic), it is of common practice to substitute
the model function :math:`h` by a polynomial approximation over its
whole domain of definition. Actually this approach is
strictly equivalent to:

#. Regarding the input parameters as random uniform random variables

#. Expanding any quantity of interest provided by the model onto a PC
   expansion made of Legendre polynomials

.. topic:: API:

    - See :class:`~openturns.FunctionalChaosAlgorithm`
    - See :class:`~openturns.HaarWaveletFactory`
    - See :class:`~openturns.FourierSeriesFactory`
    - See :class:`~openturns.SoizeGhanemFactory`
    - See :class:`~openturns.OrthogonalUniVariatePolynomialFamily`
    - See :class:`~openturns.OrthogonalUniVariatePolynomialFactory`


.. topic:: Examples:

    - See :doc:`/auto_meta_modeling/polynomial_chaos_metamodel/plot_functional_chaos`
    - See :doc:`/auto_functional_modeling/univariate_functions/plot_createUnivariateFunction`


.. topic:: References:

    - [lemaitre2010]_
    - [sullivan2015]_, chapter 11 section 11.3 page 237
    - [xiu2010]_
    - [soizeghanem2004]_

