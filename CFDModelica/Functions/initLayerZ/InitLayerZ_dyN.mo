within CFDModelica.Functions.initLayerZ;
function InitLayerZ_dyN
  "function that compute dy dimension for each Layer Z in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Real yFrac[J] "fraction of the base for each volume";
  input Modelica.SIunits.Distance width "width dimension of the room";
  output Modelica.SIunits.Distance dy[I,J,K - 1] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for k in 1:K - 1 loop
      for j in 1:J-1 loop
        dy[i, j, k] := (yFrac[j] + yFrac[j + 1])*0.5*width;
      end for;
      dy[i, J, k] := yFrac[J]*0.5*width;
    end for;
  end for;
end InitLayerZ_dyN;
