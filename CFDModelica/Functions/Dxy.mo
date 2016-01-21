within CFDModelica.Functions;
function Dxy
  input Integer I;
  input Modelica.SIunits.Distance dx[:];
  output Modelica.SIunits.Distance x[I + 1];
algorithm
  x[1]:= dx[1]/2;
  for i in 2:I loop
    x[i] := (dx[i-1]+dx[i])/2;
  end for;
  x[I+1] := dx[I]/2;
end Dxy;
