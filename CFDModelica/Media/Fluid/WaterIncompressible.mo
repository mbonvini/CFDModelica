within CFDModelica.Media.Fluid;
package WaterIncompressible "Water INCOMPRESSIBLE"
  extends CFDModelica.Media.Fluid.PartialModels.PartialIncompressibleFluid(
    cp_const=4186,
    rho0=1000,
    T_min=0,
    T_max=100,
    reference_T=273.15 + 25,
    reference_p=101325);
end WaterIncompressible;
