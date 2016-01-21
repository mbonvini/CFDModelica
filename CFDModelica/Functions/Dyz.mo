within CFDModelica.Functions;
function Dyz
  input Integer J;
  input Modelica.SIunits.Distance dy[:];
  output Modelica.SIunits.Distance y[J + 1];
algorithm
  y[1]:= dy[1]/2;
  for j in 2:J loop
    y[j] := (dy[j-1]+dy[j])/2;
  end for;
  y[J+1] := dy[J]/2;
end Dyz;
