within CFDModelica.Functions;
function arr3d1diffFlag "function that instantiate a 3d matrix with the same value except for
one element that has a different value "
  input Integer I "size on i-axis";
  input Integer J "size on j-axis";
  input Integer K "size on k-axis";
  input Boolean valall "common value";
  input Integer ione "i position for different element";
  input Integer jone "j position for different element";
  input Integer kone "k position for different element";
  input Boolean valone "different value";
  output Boolean arr[I,J,K] "output value";
algorithm
  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K loop
        arr[i, j, k] := valall;
      end for;
    end for;
  end for;
  arr[ione,jone,kone] := valone;
end arr3d1diffFlag;
