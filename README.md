# Bimodal Battery Electrode Optimization

A MATLAB-based computational modeling project investigating the
optimal particle-size distribution of bimodal battery electrodes
by balancing volumetric capacity and specific power.

## Project Overview

This project evaluates mixtures of 4 μm and 8 μm spherical
particles using MATLAB-based 3D particle packing simulations.

Eleven different particle compositions were simulated, with the
small-particle fraction varying from 0.0 to 1.0 in increments of 0.1.

The particle packing results were used to calculate effective
material density and volumetric capacity. The capacity results
were then combined with specific power through Min-Max
normalization to identify an optimal balance between the two
performance metrics.

## 3D Particle Packing Simulation

The 3D model places particles while satisfying minimum spacing
requirements based on particle size. Particle packing density was
evaluated for each small-to-large particle composition.

![3D Particle Packing Results](figures/particle_packing.png)

## Power-Capacity Optimization

The simulation results showed that maximum volumetric capacity
and maximum specific power occur at different particle
compositions.

Min-Max normalization was therefore applied to both performance
metrics to determine the composition providing the best overall
balance.

![Power-Capacity Optimization](figures/optimization_result.png)

## Key Results

| Performance Metric | Small Particle Ratio |
|---|---:|
| Maximum volumetric capacity | 0.4 |
| Maximum specific power | 1.0 |
| Optimal power-capacity balance | **0.5** |

## Methods

- 3D particle packing simulation
- Particle volume fraction calculation
- Effective material density analysis
- Volumetric capacity analysis
- Specific power analysis
- Min-Max normalization
- Multi-objective optimization
- Data visualization

## Tools

- MATLAB
- 3D Computational Modeling
- Data Analysis and Visualization
