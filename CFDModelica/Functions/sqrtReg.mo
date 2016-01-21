within CFDModelica.Functions;
function sqrtReg
  input Real x;
  input Real delta=0.001 "Range of significant deviation from sqrt(x)";
  output Real y;
algorithm
  y := x/sqrt(sqrt(x*x+delta*delta));
end sqrtReg;
