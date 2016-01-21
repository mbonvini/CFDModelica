within CFDModelica.Functions.initLayerZ;
function InitLayerZ_dxE
  "function that compute dx dimension for each Layer Z in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real xFrac[I] "fraction of the base for each volume";
  input Modelica.SIunits.Distance base "base dimension of the room";
  output Modelica.SIunits.Distance dx[I,J,K - 1] "dx dimension of each volume";
algorithm
  for j in 1:J loop
    for k in 1:K - 1 loop
      for i in 1:I-1 loop
        dx[i, j, k] := (xFrac[i] + xFrac[i + 1])*0.5*base;
      end for;
      dx[I, j, k] := xFrac[I]*0.5*base;
    end for;
  end for;
end InitLayerZ_dxE;
