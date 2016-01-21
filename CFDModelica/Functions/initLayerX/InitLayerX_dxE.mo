within CFDModelica.Functions.initLayerX;
function InitLayerX_dxE
  "function that compute dx dimension for each Layer X in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real xFrac[I] "fraction of the base for each volume";
  input Modelica.SIunits.Distance base "base dimension of the room";
  output Modelica.SIunits.Distance dx[I - 1,J,K] "dx dimension of each volume";
algorithm
  for j in 1:J loop
    for k in 1:K loop
      for i in 1:I-1 loop
        dx[i, j, k] := xFrac[i + 1]*base;
      end for;
    end for;
  end for;
end InitLayerX_dxE;
