within CFDModelica.Functions;
function wallDistanceX
  input Integer I "number of volumes along x axes";
  input Integer K "number of volumes along z axes";
  input Integer i "x coordinate of the x-volume";
  input Integer k "y coordinate of the x-volume";
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
    d_left   := 0;
    d_right  := 0;
    d_top    := 0;
    d_bottom := 0;
  else
    d_left   := Base*sum(X_frac[j] for j in 1:i-1);
    d_right  := Base*sum(X_frac[j] for j in i:I);
    if k == 2 then
      d_bottom := Height*(Z_frac[1]*0.5);
      d_top    := Height - d_bottom;
    elseif k==K+1 then
      d_top    := Height*(Z_frac[K]*0.5);
      d_bottom := Height - d_top;
    else
      d_bottom := Height*(Z_frac[1]*0.5 + sum(0.5*(Z_frac[j-1]+Z_frac[j]) for j in 2:k-1));
      d_top    := Height - d_bottom;
    end if;
  end if;
  if method == "Xu" then
    length := min(d_top,d_bottom);
  else
    length := min(min(min(d_left,d_right),d_top),d_bottom);
  end if;
end wallDistanceX;
