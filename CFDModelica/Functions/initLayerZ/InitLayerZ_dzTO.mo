within CFDModelica.Functions.initLayerZ;
function InitLayerZ_dzTO
  "function that compute dz dimension for each Layer Z in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real zFrac[K] "fraction of the base for each volume";
  input Modelica.SIunits.Distance height "width dimension of the room";
  output Modelica.SIunits.Distance dz[I,J,K - 1] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        dz[i, j, k] := zFrac[k + 1]*
          height;
      end for;
    end for;
  end for;
end InitLayerZ_dzTO;
