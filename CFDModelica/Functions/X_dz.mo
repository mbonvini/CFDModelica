within CFDModelica.Functions;
function X_dz
  input Integer I;
  input Integer J;
  input Integer K;
  input Real Z_frac[K];
  input Real height;
  output Real dz[I-1,J,K];
algorithm
  for i in 1:I-1 loop
    for j in 1:J loop
      for k in 1:K loop
        dz[i, j, k] := 5;
                       //Z_frac[k]*height;
      end for;
    end for;
  end for;
end X_dz;
