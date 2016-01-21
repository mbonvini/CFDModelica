within CFDModelica.Functions;
function wallDistanceY_3D
  input Integer I "number of volumes along x axes";
  input Integer J "number of volumes along y axes";
  input Integer K "number of volumes along z axes";
  input Integer i "x coordinate of the y-volume";
  input Integer j "y coordinate of the y-volume";
  input Integer k "y coordinate of the y-volume";
  input String method = "Xu"
    "String that indicated the measure to be implemented";
  input Real X_frac[I];
  input Real Y_frac[J];
  input Real Z_frac[K];
  input Real Base;
  input Real Width;
  input Real Height;
  output Real length;
protected
  Real d_left;
  Real d_right;
  Real d_front;
  Real d_rear;
  Real d_top;
  Real d_bottom;
algorithm
  if (i==1 or i==I+2 or j==1 or j==J+2 or k==1 or k==K+2) then
    d_left   := 0;
    d_right  := 0;
    d_front  := 0;
    d_rear   := 0;
    d_top    := 0;
    d_bottom := 0;
  else

    // distance along x direction
    if i == 2 then
      d_left  := Base*(X_frac[1]*0.5);
      d_right := Base - d_left;
    elseif i==I+1 then
      d_right := Base*(X_frac[I]*0.5);
      d_left  := Base - d_right;
    else
      d_left  := Base*(X_frac[1]*0.5 + sum(0.5*(X_frac[l-1]+X_frac[l]) for l in 2:i-1));
      d_right := Base - d_left;
    end if;

    // distance along y direction
    d_front := Width*sum(Y_frac[l] for l in 1:j-1);
    d_rear  := Width*sum(Y_frac[l] for l in j:J);

    // distance along z direction
    if k == 2 then
      d_bottom := Height*(Z_frac[1]*0.5);
      d_top    := Height - d_bottom;
    elseif k==K+1 then
      d_top    := Height*(Z_frac[K]*0.5);
      d_bottom := Height - d_top;
    else
      d_bottom := Height*(Z_frac[1]*0.5 + sum(0.5*(Z_frac[l-1]+Z_frac[l]) for l in 2:k-1));
      d_top    := Height - d_bottom;
    end if;

  end if;

  // the comparison depends on the method
  if method == "Xu" then
    length := min([d_top,d_bottom,d_left,d_right]);
  else
    length := min([d_left,d_right,d_top,d_bottom,d_front,d_rear]);
  end if;

end wallDistanceY_3D;
