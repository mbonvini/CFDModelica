within CFDModelica.Functions;
function MU_walls
  input Modelica.SIunits.DynamicViscosity mu
    "dynamic viscosity for the viscous case";
  input Real Kturb "turbulent constant";
  input Modelica.SIunits.Velocity V "air velocity";
  input Modelica.SIunits.Density rho "air reference density";
  input Modelica.SIunits.Distance y_wall "distance from the wall";
  output Modelica.SIunits.DynamicViscosity mu_wall
    "dynamic viscosity of air in the near wall zone";
  parameter Real Y_PLUS_min = 5
    "minimum value for the normalised wall distance";
  parameter Real Y_PLUS_max = 5
    "minimum value for the normalised wall distance";
  parameter Real E = 8.9 "log-law model constant";
  parameter Real k = 0.41 "Von-Karman constant";
protected
  Real mu_t "viscosity in the turbulent case";
  Real mu_Y_PLUS_max "value for the end of the linear phase";
  Real y_plus "normalised distance";
algorithm
  // compute the wall viscosity
  mu_t := Kturb*sqrtReg(V^2, 1e-8)*rho*y_wall;

  // compute y_plus
  y_plus := mu_t/mu/k;

  // maximum value of mu after the linear phase
  mu_Y_PLUS_max := mu_t/Modelica.Math.log(E*Y_PLUS_max);

  // wall viscosity
  if y_plus < Y_PLUS_min then
      mu_wall := mu;
  elseif y_plus >= Y_PLUS_min and y_plus <= Y_PLUS_max then
      mu_wall := mu + (y_plus - Y_PLUS_min)*(mu_Y_PLUS_max - mu)/(Y_PLUS_max - Y_PLUS_min);
  else
      mu_wall := mu_t/Modelica.Math.log(E*y_plus);
  end if;

end MU_walls;
