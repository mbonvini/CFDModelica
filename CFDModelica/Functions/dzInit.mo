within CFDModelica.Functions;
function dzInit "function that compute dx dimension for each volume in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real zFrac[K] "fraction of the base for each volume";
  input Modelica.SIunits.Distance height "height dimension of the room";
  output Modelica.SIunits.Distance dz[I,J,K] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K loop
        dz[i, j, k] := zFrac[k]*height;
      end for;
    end for;
  end for;
end dzInit;
