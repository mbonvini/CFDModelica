within CFDModelica.Media.Fluid.PartialModels;
partial package PartialSimpleLinearizedIdealGas "Medium model of Ideal gas with constant cp and cv. 
    All other quantities, e.g. transport properties, are constant."
  extends Modelica.Media.Interfaces.PartialPureSubstance(final singleState = false);
  import SI = Modelica.SIunits;
  constant SI.SpecificHeatCapacity cp_const = cv_const + R_gas
  "Constant specific heat capacity at constant pressure";
  constant SI.SpecificHeatCapacity cv_const= 1005.45
  "Constant specific heat capacity at constant volume";
  constant SI.SpecificHeatCapacity R_gas "medium specific gas constant";
  constant SI.MolarMass MM_const "Molar mass";
  constant SI.DynamicViscosity eta_const "Constant dynamic viscosity";
  constant SI.ThermalConductivity lambda_const "Constant thermal conductivity";
  constant SI.Temperature T_min "Minimum temperature valid for medium model";
  constant SI.Temperature T_max "Maximum temperature valid for medium model";
  constant SI.Temperature T0=0
  "Zero enthalpy temperature, temperature linearization value";
  constant SI.Temperature T0lin = reference_T;
  constant SI.Pressure p0 = reference_p "Pressure linearization value";
  constant SI.Pressure p_ref=reference_p "Absolute pressure reference";
  constant SI.Density d0 "Density linearization value";
  constant Real ComprCoeff "Air comprimibility coefficient";
  constant Real ThermalExpCoeff "Air thermal expansion coefficient";


  redeclare record extends ThermodynamicState "thermodynamic state"
    SI.Pressure p "Relative pressure of medium";
    SI.Temperature T "Temperature of medium";
  end ThermodynamicState;


  replaceable function specificInternalEnergy_ph
  "Return specific internal energy from p and h"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Pressure p "Pressure";
    input Modelica.SIunits.SpecificEnthalpy h "Specific enthalpy";
    output Modelica.SIunits.SpecificInternalEnergy u "Specific internal energy";
  algorithm
    u := cv_const*(h/cp_const + T0);
  end specificInternalEnergy_ph;


  replaceable function specificInternalEnergy_pT
  "Return specific internal energy from p and T"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.Pressure p "Pressure";
    input SI.Temperature T "Temperature";
    output Modelica.SIunits.SpecificInternalEnergy u "specific internal energy";
  algorithm
    u := cv_const*T;
  end specificInternalEnergy_pT;


  redeclare function temperature_phX
  "Return temperature from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input SI.Pressure p "Pressure";
    input SI.SpecificEnthalpy h "Specific enthalpy";
    input SI.MassFraction X[nX] "Mass fractions";
    output SI.Temperature T "Temperature";
  algorithm
    T := h/cp_const + T0;
  end temperature_phX;


  redeclare function temperature_ph "Return temperature from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input SI.Pressure p "Pressure";
    input SI.SpecificEnthalpy h "Specific enthalpy";
    output SI.Temperature T "Temperature";
  algorithm
    T := h/cp_const + T0;
  end temperature_ph;


  redeclare function density_phX "Return density from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input SI.Pressure p "Pressure";
    input SI.SpecificEnthalpy h "Specific enthalpy";
    input SI.MassFraction X[nX] "Mass fractions";
    output SI.Density d "density";
  algorithm
    d := density(setState_phX(p,h,X));
  end density_phX;


  redeclare function setState_pTX
  "Return thermodynamic state from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input SI.Pressure p "Pressure";
    input SI.Temperature T "Temperature";
    input SI.MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state := ThermodynamicState(p=p,T=T);
  end setState_pTX;


  redeclare function setState_phX
  "Return thermodynamic state from p, h, and X or Xi"
    extends Modelica.Icons.Function;
    input SI.Pressure p "Pressure";
    input SI.SpecificEnthalpy h "Specific enthalpy";
    input SI.MassFraction X[:]=reference_X "Mass fractions";
    output ThermodynamicState state;
  algorithm
    state := ThermodynamicState(p=p,T=T0+h/cp_const);
  end setState_phX;


  redeclare function extends pressure "return pressure of ideal gas"
  algorithm
    p := state.p;
  end pressure;


  redeclare function extends temperature "return temperature of ideal gas"
  algorithm
    T := state.T;
  end temperature;


  redeclare function extends density "return density of ideal gas"
  algorithm
   d := d0 + ComprCoeff*(state.p) - ThermalExpCoeff*(state.T-T0lin);
  end density;


  redeclare function extends specificEnthalpy "Return specific enthalpy"
      extends Modelica.Icons.Function;
  algorithm
    Modelica.Constants.h := cp_const*(state.T - T0);
  end specificEnthalpy;


  redeclare function extends specificInternalEnergy
  "Return specific internal energy"
    extends Modelica.Icons.Function;
  algorithm
    u := cv_const*state.T;
  end specificInternalEnergy;


  redeclare function extends specificEntropy "Return specific entropy"
      extends Modelica.Icons.Function;
  algorithm
    s := cp_const*Modelica.Math.log(state.T/T0) - R_gas*Modelica.Math.log((state.p + p_ref)/reference_p);
  end specificEntropy;


  redeclare function extends dynamicViscosity "Return const dynamic viscosity"
  algorithm
    eta := eta_const;
  end dynamicViscosity;


  redeclare function extends thermalConductivity "Return thermal conductivity"
  algorithm
    lambda := lambda_const;
  end thermalConductivity;


  redeclare function extends specificHeatCapacityCp
  "Return specific heat capacity at constant pressure"
  algorithm
    cp := cp_const;
  end specificHeatCapacityCp;


  redeclare function extends specificHeatCapacityCv
  "Return specific heat capacity at constant volume"
  algorithm
    cv := cv_const;
  end specificHeatCapacityCv;


  redeclare function extends isentropicExponent "Return isentropic exponent"
  algorithm
    Modelica.Constants.gamma := cp_const/cv_const;
  end isentropicExponent;


  redeclare function extends velocityOfSound "Velocity of sound "
  algorithm
    a := sqrt(cp_const/cv_const*R_gas*state.T);
  end velocityOfSound;


  redeclare function specificEnthalpy_pTX
  "Return specific enthalpy from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input SI.Pressure p "Pressure";
    input SI.Temperature T "Temperature";
    input SI.MassFraction X[nX] "Mass fractions";
    output SI.SpecificEnthalpy h "Specific enthalpy at p, T, X";
  algorithm
    h := cp_const*(T-T0);
  end specificEnthalpy_pTX;

end PartialSimpleLinearizedIdealGas;
