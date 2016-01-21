within CFDModelica.Functions;
function Dzy
  input Integer K;
  input Modelica.SIunits.Distance dz[:];
  output Modelica.SIunits.Distance z[K + 1];
algorithm
  z[1]:= dz[1]/2;
  for k in 2:K loop
    z[k] := (dz[k - 1] + dz[k])
      /2;
  end for;
  z[K+1] := dz[K]/2;
end Dzy;
