within CFDModelica.Functions;
function obstacleY
  "function that compute dx dimension for each volume in a room"
  input Integer I "number of volumes along x-axis";
  input Integer J "number of volumes along y-axis";
  input Integer K "number of volumes along z-axis";
  input Integer Io;
  input Integer Ko;
  output Boolean fixed[I,J-1,K] "dx dimension of each volume";
algorithm
  for i in 1:I loop
    for j in 1:J-1 loop
      for k in 1:K loop
        if k <= Ko and i == Io then
          fixed[i, j, k] := true;
        else
          fixed[i, j, k] := false;
        end if;
      end for;
    end for;
  end for;
end obstacleY;
