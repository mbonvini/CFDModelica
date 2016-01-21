within CFDModelica.Functions;
function H_walls
  input Modelica.SIunits.ThermalConductivity lambda
    "thermal conductivity for the viscous case";
  input Modelica.SIunits.DynamicViscosity mu
    "dynamic viscosity for the viscous case";
  input Modelica.SIunits.DynamicViscosity mu_t
    "dynamic viscosity for the turbulent case";
  input Modelica.SIunits.Length d "reference distance";
  input Real h "convective heat trasfer coeff.";
  input Boolean fixed
    "flag that indicates a fixed convective heat trasfer coeff.";
  output Modelica.SIunits.ThermalConductivity lambda_wall
    "wall thermal conductivity due to the presence of air flow";
  parameter Real Y_PLUS_min = 5
    "minimum value for the normalised wall distance";
  parameter Real Y_PLUS_max = 5
    "minimum value for the normalised wall distance";
  parameter Real E = 8.9 "log-law model constant";
  parameter Real k = 0.41 "Von-Karman constant";
protected
  Real Lambda_Y_PLUS_max "value for the end of the linear phase";
  Real y_plus "normalised distance";
algorithm

  if fixed then
    // fixed h
    lambda_wall := h*d;
  else
    // compute the heat transfer coefficient
      y_plus := mu_t/mu/k;

      Lambda_Y_PLUS_max := lambda*Y_PLUS_max*k/(Modelica.Math.log(E*Y_PLUS_max));

      if y_plus < Y_PLUS_min then
        lambda_wall := lambda;
      elseif y_plus >= Y_PLUS_min and y_plus <= Y_PLUS_max then
        lambda_wall := lambda + (y_plus - Y_PLUS_min)*(Lambda_Y_PLUS_max - lambda)/(Y_PLUS_max - Y_PLUS_min);
      else
        lambda_wall := lambda*y_plus*k/(Modelica.Math.log(E*y_plus));
      end if;

  end if;

end H_walls;
