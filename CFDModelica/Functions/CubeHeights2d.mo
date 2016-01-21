within CFDModelica.Functions;
function CubeHeights2d "function that instantiate a 3d matrix with the same value except for
one element that has a different value "
  input Integer I "size on i-axis";
  input Integer K "size on k-axis";
  input Modelica.SIunits.Length dz[K] "relative height of each volume";
  output Modelica.SIunits.Length h[I,1,K] "output value";
protected
  Modelica.SIunits.Length H[K];
  Modelica.SIunits.Length Htot;
algorithm
  H[1] := dz[1]/2;
  Htot := dz[1];
  for k in 2:K loop
    H[k] := dz[k]/2 + Htot;
    Htot := Htot + dz[k];
  end for;
  for k in 1:K loop
    for i in 1:I loop
      h[i, 1, k] := Htot - H[k];
    end for;
  end for;
end CubeHeights2d;
