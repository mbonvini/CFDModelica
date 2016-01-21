within CFDModelica.Functions.initLayerY;
function InitLayerY_dzBO
  "function that compute dz dimension for each Layer Y in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real zFrac[K] "fraction of the base for each volume";
  input Modelica.SIunits.Distance height "width dimension of the room";
  output Modelica.SIunits.Distance dz[I,J - 1,K] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for j in 1:J-1 loop
      dz[i,j,1]:= zFrac[1]*0.5*height;
      for k in 2:K loop
        dz[i, j, k] := (zFrac[k - 1] +
          zFrac[k])*0.5*height;
      end for;
    end for;
  end for;
end InitLayerY_dzBO;
