within CFDModelica.Functions.initLayerX;
function InitLayerX_dzBO
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
      dz[i,j,1]:= zFrac[1]*0.5*height;
      for k in 2:K loop
        dz[i, j, k] := (zFrac[k - 1] +
          zFrac[k])*0.5*height;
      end for;
    end for;
  end for;
end InitLayerX_dzBO;
