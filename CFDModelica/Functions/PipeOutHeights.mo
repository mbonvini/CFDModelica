within CFDModelica.Functions;
function PipeOutHeights
  input Integer k_pipes;
  input Real dz_out[k_pipes];
  input Modelica.SIunits.Length offset "offset";
  output Real h_out[k_pipes];
algorithm
  h_out[1] := dz_out[1]/2;

  for i in 2:k_pipes loop
    h_out[i] := h_out[i-1] + dz_out[i-1]*0.5 + dz_out[i]*0.5;
  end for;

end PipeOutHeights;
