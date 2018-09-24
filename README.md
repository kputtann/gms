# gms
A Generalized H-infinity Mixed Sensitivity Convex Approach to Multivariable Control Design Subject to Simultaneous Output and Input Loop-Breaking Specifications

Computes a H-infinity based Feedback Controller based on multiobjective constrained convex optimization.

Outline of steps for GMS problem setup:
   - Form the design plant:
       - Define the original plant
       - Integrator augmentation if needed
       - Bilinear transformation values if needed

   - Select weighting functions:
       - Tradeoff param rho
       - W for obj
       - W for constraint

   - Select optimization params:
       - LB and UB
       - Init point
       - Maximum number of iterations

   - Select Youla/Zames parametrization:
       - Select Youla or Zames
       - Initial controller

   - Finite Dimensionality
       - Basis params

   - Objective function:
           - sum/max/stacking

   - Find initial controller (Ko, F, L)

   - Youla parameterization

   - Find Initial Q parameter using initial controller (Ko, F, L)

   - Extract required data from problem setup

   - Vectorize the optimization problem

   - Optimization process
       - define how subgradient is picked based on sum/max/stacking

   - form Q using the optimized variables and bases

   - form Controller K using the obtained Q

   - Inverse bilinear transformation if needed

   - Inverse of integrator augmentation if needed

   - Compute OL and CL maps
