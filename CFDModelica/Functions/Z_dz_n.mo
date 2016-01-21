within CFDModelica.Functions;
function Z_dz_n
  input Integer I;
  input Integer J;
  input Integer K;
  input Real Z_frac[K];
  input Real height;
  output Real dz_n[I,J,K-1];
algorithm
  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        dz_n[i, j, k] := height*Z_frac[k
           + 1];
      end for;
    end for;
  end for;
end Z_dz_n;
