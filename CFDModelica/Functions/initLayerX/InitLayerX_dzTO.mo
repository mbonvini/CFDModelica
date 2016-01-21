within CFDModelica.Functions.initLayerX;
function InitLayerX_dzTO
  "function that compute dz dimension for each Layer X in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real zFrac[K] "fraction of the base for each volume";
  input Modelica.SIunits.Distance height "width dimension of the room";
  output Modelica.SIunits.Distance dz[I - 1,J,K] "dx dimension of each volume";
algorithm
  for i in 1:I-1 loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        dz[i, j, k] := (zFrac[k] +
          zFrac[k + 1])*0.5*height;
      end for;
      dz[i,j,K]:= zFrac[K]*0.5*height;
    end for;
  end for;
end InitLayerX_dzTO;
