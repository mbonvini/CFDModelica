within CFDModelica.Functions;
function PipeInHeights
  input Integer K;
  input Integer k_pipes;
  input Real roomHeight;
  input Real dz_in[k_pipes];
  input Modelica.SIunits.Length offset "offset";
  output Real h_in[k_pipes];
algorithm
  h_in[1] := roomHeight + dz_in[1]/2 + offset;

  for i in 2:k_pipes loop
    h_in[i] := h_in[i-1] + dz_in[i-1]*0.5 + dz_in[i]*0.5;
  end for;
end PipeInHeights;
