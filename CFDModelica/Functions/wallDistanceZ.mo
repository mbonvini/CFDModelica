within CFDModelica.Functions;
function wallDistanceZ
  input Integer I "number of volumes along x axes";
  input Integer K "number of volumes along z axes";
  input Integer i "x coordinate of the z-volume";
  input Integer k "y coordinate of the z-volume";
  input String method = "Xu"
    "String that indicated the measure to be implemented";
  input Real X_frac[I];
  input Real Z_frac[K];
  input Real Base;
  input Real Height;
  output Real length;
protected
    Real d_left;
    Real d_right;
    Real d_top;
    Real d_bottom;
algorithm
  if (i==1 or i==I+2 or k==1 or k==K+2) then
    d_left   :=0;
    d_right  := 0;
    d_top    := 0;
    d_bottom := 0;
  else
    d_top    := Height*sum(Z_frac[j] for j in k:K);
    d_bottom := Height*sum(Z_frac[j] for j in 1:k-1);
    if i == 2 then
      d_left  := Base*(X_frac[1]*0.5);
      d_right := Base - d_left;
    elseif i==I+1 then
      d_right := Base*(X_frac[I]*0.5);
      d_left  := Base - d_right;
    else
      d_left  := Base*(X_frac[1]*0.5 + sum(0.5*(X_frac[j-1]+X_frac[j]) for j in 2:i-1));
      d_right := Base - d_left;
    end if;
  end if;
  if method == "Xu" then
    length := min(d_left,d_right);
  else
    length := min(min(min(d_left,d_right),d_top),d_bottom);
  end if;
end wallDistanceZ;
