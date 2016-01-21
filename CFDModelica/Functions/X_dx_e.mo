within CFDModelica.Functions;
function X_dx_e
  input Integer I;
  input Integer J;
  input Integer K;
  input Real X_frac[I];
  input Real base;
  output Real dx_e[I-1,J,K];
algorithm
  for i in 1:I-1 loop
    for j in 1:J loop
      for k in 1:K loop
        dx_e[i, j, k] := X_frac[i + 1]*base;
      end for;
    end for;
  end for;
end X_dx_e;
