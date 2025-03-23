.. barotropy documentation master file, created by
   sphinx-quickstart on Thu Sep 28 17:55:22 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to barotropy's documentation!
=====================================

.. module:: barotropy
   :synopsis: A Python package to generate barotropic fluid property models.

``ThermOpt`` is a Python package for the design and analysis of thermodynamic cycles.
 generate barotropic fluid property models. 

The thermodynamic properties are computed using the CoolProp library.

The fluid properties between the saturation line and the spinodal line can be 
computed according to phase-equilibrium or by extrapolating the Helmholtz-energy 
equation of state into the metastable region. Check the 
`documentation <./documentation/thermodynamic_properties.rst>`_ for more information 
about extrapolating the equation of state beyond the saturation line and its limitations.

Use the panel to the left or the table of contents below to navigate the documentation.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   source/installation
   source/tutorials
   source/cycle_optimization
   source/developer_guide
   source/bibliography
   source/api/thermopt

