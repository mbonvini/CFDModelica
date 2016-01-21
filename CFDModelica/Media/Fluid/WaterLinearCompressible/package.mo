within CFDModelica.Media.Fluid;
package WaterLinearCompressible "Water: Simple linear compressible model (0..80 degC)"
  extends CFDModelica.Media.Fluid.PartialModels.PartialSimpleLinearisedFluid(
  cp_const=4184,
  rho0=997.047,
  T_min=Cv.from_degC(0),
  T_max=Cv.from_degC(100),
  reference_T=273.15 + 25,
  reference_p=101325,
  eta_const=1e-3,
  lambda_const=0.598,
  ComprCoeff=2e-6,
  ThermalExpCoeff=0.20907306);
  import Cv = Modelica.SIunits.Conversions;
  import Modelica.Constants;
end WaterLinearCompressible;
