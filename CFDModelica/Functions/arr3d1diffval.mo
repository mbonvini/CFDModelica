within CFDModelica.Functions;
function arr3d1diffval "function that instantiate a 3d matrix with the same value except for
one element that has a different value "
  input Integer I "size on i-axis";
  input Integer J "size on j-axis";
  input Integer K "size on k-axis";
  input Real valall "common value";
  input Integer ione "i position for different element";
  input Integer jone "j position for different element";
  input Integer kone "k position for different element";
  input Real valone "different value";
  output Real arr[I,J,K] "output value";
algorithm
  arr := valall*ones(I,J,K);
  arr[ione,jone,kone] := valone;
end arr3d1diffval;
