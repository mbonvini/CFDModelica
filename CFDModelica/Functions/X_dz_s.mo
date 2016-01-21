within CFDModelica.Functions;
function X_dz_s
  input Integer I;
  input Integer J;
  input Integer K;
  input Real Z_frac[K];
  input Real height;
  output Real dz_s[I-1,J,K];
algorithm
  for i in 1:I-1 loop
    for j in 1:J loop
      dz_s[i,j,1] := Z_frac[1]*height/2;
    end for;
  end for;
  for i in 1:I-1 loop
    for j in 1:J loop
      for k in 2:K loop
        dz_s[i, j, k] := (Z_frac[k] +
          Z_frac[k - 1])*height/2;
      end for;
    end for;
  end for;
end X_dz_s;
