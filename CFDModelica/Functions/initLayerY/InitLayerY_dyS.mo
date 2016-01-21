within CFDModelica.Functions.initLayerY;
function InitLayerY_dyS
  "function that compute dy dimension for each Layer Y in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real yFrac[J] "fraction of the base for each volume";
  input Modelica.SIunits.Distance width "width dimension of the room";
  output Modelica.SIunits.Distance dy[I,J - 1,K] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for k in 1:K loop
      for j in 1:J-1 loop
        dy[i, j, k] := yFrac[j]*width;
      end for;
    end for;
  end for;
end InitLayerY_dyS;
