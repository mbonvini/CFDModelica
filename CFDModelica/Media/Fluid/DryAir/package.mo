within CFDModelica.Media.Fluid;
package DryAir "Air: Simple dry air model (0..100 degC)"
  extends CFDModelica.Media.Fluid.PartialModels.PartialSimpleLinearizedIdealGas(
  mediumName="SimpleAir",
  cv_const=1005.45,
  MM_const=0.0289651159,
  R_gas=Constants.R/0.0289651159,
  eta_const=1.82e-5,
  lambda_const=0.026,
  T_min=Cv.from_degC(0),
  T_max=Cv.from_degC(100),
  ComprCoeff=1/(R_gas*T0lin),
  ThermalExpCoeff=p0/(R_gas*T0lin^2),
  d0=p0/(R_gas*T0lin),
  reference_T=273.15 + 15,
  reference_p=101321);
  import Cv = Modelica.SIunits.Conversions;
  import Modelica.Constants;

  annotation (Documentation(info="<html>
                              <h2>Simple Ideal gas air model for low temperatures<h1>
                              <p>This model demonstrats how to use the PartialSimpleIdealGas base class to build a
                              simple ideal gas model with a limited temperature validity range.</p>
                              </html>"));
end DryAir;
