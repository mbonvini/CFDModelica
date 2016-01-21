within CFDModelica.Functions;
function Z_dx_w
  input Integer I;
  input Integer J;
  input Integer K;
  input Real X_frac[I];
  input Real base;
  output Real dx_w[I,J,K-1];
algorithm
  for k in 1:K - 1 loop
    for j in 1:J loop
      dx_w[1, j, k] := X_frac[1]*base/2;
    end for;
  end for;
  for i in 2:I loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        dx_w[i, j, k] := (X_frac[i - 1] + X_frac[i])*base/2;
      end for;
    end for;
  end for;
end Z_dx_w;
