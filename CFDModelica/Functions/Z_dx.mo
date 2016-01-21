within CFDModelica.Functions;
function Z_dx
  input Integer I;
  input Integer J;
  input Integer K;
  input Real X_frac[I];
  input Real base;
  output Real dx[I,J,K-1];
algorithm
  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        dx[i, j, k] := X_frac[i]*base;
      end for;
    end for;
  end for;
end Z_dx;
