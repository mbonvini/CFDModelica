within CFDModelica.Functions;
function dyInit "function that compute dx dimension for each volume in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real yFrac[J] "fraction of the base for each volume";
  input Modelica.SIunits.Distance width "width dimension of the room";
  output Modelica.SIunits.Distance dy[I,J,K] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for k in 1:K loop
      for j in 1:J loop
        dy[i, j, k] := yFrac[j]*width;
      end for;
    end for;
  end for;
end dyInit;
