within CFDModelica.Media.Test;
model ProvaModelicaWater_h0
  "test - this models had been used to test the value of h0"
  package Medium = Modelica.Media.Water.StandardWater;
  parameter Modelica.SIunits.Temperature T0=Medium.reference_T;
  parameter Modelica.SIunits.Pressure p0=Medium.reference_p;
  parameter Modelica.SIunits.SpecificEnthalpy h0=
      Modelica.Media.Water.IF97_Utilities.h_pT(p0, T0);
end ProvaModelicaWater_h0;
