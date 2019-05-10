# PRRMD
Parameter Range Reduction using Monotonic Discretizations

Copyright 2019 Andrew Skelton and Allan R. Willms.

This C software source code is distributed under the GNU General Public License, Version 3.

Currently the software here is not fully tested, but should be operational.  A usage manual is yet to be written.

Bug reports and comments should be sent to Allan Willms.

<h3>Description</h3>

Consider a Differential-Algebraic Equation of the form<br>
x' = f(t,x,lambda)<br>
0 = g(t,x,lambda)<br>
y = phi(t,x,lambda)<br>
Here x are the state variables, y are the observed quantities, and lambda are the parameters.  Given a data set {(t<sub>i</sub>, y<sub>i</sub>)} and initial ranges for all of the paramters lambda, the goal is to reduce these ranges as much as possible by eliminating regions of parameter space that are inconsistent with the data.
   
PRRMD reduces initial parameter ranges using monotonic discretizations.  

The PRRMD algorithm is described in
<ul> 
  <li> Allan R. Willms, 2007. Parameter Range Reduction for ODE Models Using Cumulative Backward Differentiation Formulas, J. Comput. Appl. Math. 203 (1), pp. 87-102. 
    <li> Allan R. Willms, Emily K. Szusz, 2013. Parameter Range Reduction for ODE Models Using Monotonic Discretizations, J. Comput. Appl. Math. 247, pp. 124-151. 
      <li> Andrew Skelton, Allan R. Willms, 2015. Parameter Range Reduction in Ordinary Differential Equation Models, Journal of Scientific Computing 62 (2), pp. 517-531, doi=10.1007/s10915-014-9865-6.
 </ul>
  
