within CFDModelica.Functions;
function X_dz_n
  input Integer I;
  input Integer J;
  input Integer K;
  input Real Z_frac[K];
  input Real height;
  output Real dz_n[I-1,J,K];
algorithm
  for i in 1:I-1 loop
    for j in 1:J loop
      dz_n[i,j,K] := Z_frac[K]*height/2;
    end for;
  end for;
  for i in 1:I-1 loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        dz_n[i, j, k] := (Z_frac[k] +
          Z_frac[k + 1])*height/2;
      end for;
    end for;
  end for;
end X_dz_n;
