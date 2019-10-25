=========
Use-cases
=========

Introduction
------------

This document presents an index of the physical use-cases and 
the application of various methods on it. 

Each use-case is made of:

- a function g which computes outputs depending on inputs,
- a probabilistic model of the inputs,
- a quantity of interest.

The quantity of interest may be for example the mean, the probability of 
exceeding a threshold (reliability), Sobol' indices (sensitivity 
analysis), etc...

Cantilever beam
---------------

- :doc:`/examples/reliability_sensitivity/central_tendency`
- :doc:`/examples/reliability_sensitivity/estimate_probability_form`
- :doc:`/examples/reliability_sensitivity/estimate_probability_importance_sampling`
- :doc:`/examples/reliability_sensitivity/estimate_probability_randomized_qmc`
- :doc:`/examples/reliability_sensitivity/estimate_probability_directional_sampling`
- :doc:`/examples/meta_modeling/kriging_cantilever_beam`
- :doc:`/examples/meta_modeling/kriging_beam_trend`
- :doc:`/examples/meta_modeling/kriging_beam_arbitrary_trend`
- :doc:`/examples/meta_modeling/kriging_hyperparameters_optimization`
- :doc:`/examples/meta_modeling/chaos_cantilever_beam_integration`
- :doc:`/examples/numerical_methods/minmax_optimization`
- :doc:`/examples/numerical_methods/control_termination`
- :doc:`/examples/numerical_methods/minmax_by_random_design`

Ishigami function
-----------------

- :doc:`/examples/data_analysis/sample_correlation`
- :doc:`/examples/meta_modeling/chaos_ishigami_grouped_indices`
- :doc:`/examples/meta_modeling/chaos_ishigami`
- :doc:`/examples/reliability_sensitivity/sensitivity_sobol`

Flooding model
--------------

- :doc:`/examples/reliability_sensitivity/flood_model`
- :doc:`/examples/data_analysis/compare_unconditional_conditional_histograms`
- :doc:`/examples/calibration/calibration_flooding`
- :doc:`/examples/calibration/bayesian_calibration_flooding`

Axial stressed beam
-------------------

- :doc:`reliability_sensitivity/axial_stressed_beam_quickstart`
- :doc:`reliability_sensitivity/axial_stressed_beam`

Viscous free fall
-----------------

- :doc:`/examples/meta_modeling/viscous_fall_metamodel`
- :doc:`/examples/functional_modeling/viscous_fall_field_function`

Logistic growth
---------------

- :doc:`/examples/functional_modeling/logistic_growth_model`
- :doc:`/examples/calibration/calibration_logistic`

Deflection of a tube
--------------------

- :doc:`/examples/calibration/calibration_deflection_tube`

Chaboche mechanical model
-------------------------

- :doc:`/examples/calibration/calibration_chaboche`

Borehole
--------

- :doc:`/examples/reliability_sensitivity/functional_chaos_sensitivity`

Branin function
---------------

- :doc:`/examples/numerical_methods/ego`

Ackley function
---------------

- :doc:`/examples/numerical_methods/ego`

