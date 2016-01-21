within CFDModelica.Media.Test;
model ProvaBEPIWater_h0
  "test - this models had been used to test the value of h0"
  package Medium = CFDModelica.Media.Fluid.WaterLinearCompressible;
  parameter Modelica.SIunits.Temperature T0=Medium.reference_T;
  parameter Modelica.SIunits.Pressure p0=Medium.reference_p;
  parameter Modelica.SIunits.SpecificEnthalpy h0=
      Medium.specificEnthalpy_pT(p0, T0);
  parameter Modelica.SIunits.Density d0=Medium.density_pT(p0, T0);
end ProvaBEPIWater_h0;
