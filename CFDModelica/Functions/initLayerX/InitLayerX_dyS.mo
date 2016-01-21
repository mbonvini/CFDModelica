within CFDModelica.Functions.initLayerX;
function InitLayerX_dyS
  "function that compute dy dimension for each Layer X in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real yFrac[J] "fraction of the base for each volume";
  input Modelica.SIunits.Distance width "width dimension of the room";
  output Modelica.SIunits.Distance dy[I - 1,J,K] "dx dimension of each volume";
algorithm
  for i in 1:I-1 loop
    for k in 1:K loop
      dy[i, 1, k] := yFrac[1]*width*0.5;
      for j in 2:J loop
        dy[i, j, k] := (yFrac[j - 1] + yFrac[j])*0.5*width;
      end for;
    end for;
  end for;
end InitLayerX_dyS;
