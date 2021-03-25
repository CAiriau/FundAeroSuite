=======================
Compressible Aero Tool
=======================

Introduction
-------------

* To test the tool just run the function :mod:`run_cat` in the **CompressibleFlow** directory.

* As written before many functions can be used independantly and integrated into a new tool. 

* The functions are listed by themes.

.. automodule:: CompressibleFlow

.. autofunction:: run_cat

.. automodule:: CompressibleFlow.cat

Generalities
-------------

.. autofunction:: air_properties
.. autofunction:: compressible_flow
.. autofunction:: display_menu
.. autofunction:: sound_celerity

Isentropic flow
---------------

.. autofunction:: isentropic_evolution
.. autofunction:: isentropic_ratios
.. autofunction:: Prandtl_Meyer
.. autofunction:: inverse_Prandtl_Meyer
.. autofunction:: table_omega
.. autofunction:: epicycloide

shocks
-------

.. autofunction:: Rho2overRho1
.. autofunction:: shock_angle
.. autofunction:: solve_shock
.. autofunction:: downstream_normal_Mach
.. autofunction:: for_iso_sigma
.. autofunction:: iso_mach
.. autofunction:: oblique_shock_curve
.. autofunction:: P2overP1
.. autofunction:: Pi2overPi1
.. autofunction:: theta_max_curve
.. autofunction:: unitary_downstream_Mach_curve

Fanno flow
----------

.. autofunction:: Fanno_flow
.. autofunction:: Fanno
.. autofunction:: inverse_Fanno
.. autofunction:: Exercise_Fanno
.. autofunction:: Fanno_entropy
.. autofunction:: Mass_Flow_Rate
.. autofunction:: qm_Mach

Rayleigh flow
-------------

.. autofunction:: inverse_Rayleigh
.. autofunction:: plot_Rayleigh_entropy
.. autofunction:: Rayleigh_flow
.. autofunction:: Rayleigh

Misc
-----

.. autofunction:: hline













