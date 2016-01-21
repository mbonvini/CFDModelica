within CFDModelica.Flow3D;
partial model Room3D_CFDwall_base
  /*********** MEDIUM MODEL ************************************************************************/
  replaceable package Medium = CFDModelica.Media.Fluid.DryAir
    "|Media| Medium model";
  /*********** INLET AND OUTLETS *******************************************************************/
  parameter Boolean all_left = false
    "This flag enables all left fluid connectors"
  annotation (Evaluate = true,Dialog(tab="Inlet/outlet", group="Left"),choices(__Dymola_checkBox=true));
  parameter Integer jk_left[:,2] = [0,0]
    "coordinates of the (left) fluid connectors [j1,k1;j2,k2;...]"
  annotation ( Evaluate = true,Dialog(enable = not all_left,tab="Inlet/outlet", group="Left"));
  parameter Boolean all_right = false
    "This flag enables all right fluid connectors"
  annotation (Evaluate = true,Dialog(tab="Inlet/outlet", group="Right"),choices(__Dymola_checkBox=true));
  parameter Integer jk_right[:,2] = [0,0]
    "coordinates of the (right) fluid connectors [j1,k1;j2,k2;...]"
  annotation (Evaluate = true,Dialog(enable = not all_right,tab="Inlet/outlet", group="Right"));
  parameter Boolean all_front = false
    "This flag enables all front fluid connectors"
  annotation (Evaluate = true,Dialog(tab="Inlet/outlet", group="Front"),choices(__Dymola_checkBox=true));
  parameter Integer ik_front[:,2] = [0,0]
    "coordinates of the (front) fluid connectors [i1,k1;i2,k2;...]"
  annotation (Evaluate = true,Dialog(enable = not all_front,tab="Inlet/outlet", group="Front"));
  parameter Boolean all_rear = false
    "This flag enables all rear fluid connectors"
  annotation (Evaluate = true,Dialog(tab="Inlet/outlet", group="Rear"),choices(__Dymola_checkBox=true));
  parameter Integer ik_rear[:,2] = [0,0]
    "coordinates of the (rear) fluid connectors [i1,k1;i2,k2;...]"
  annotation (Evaluate = true,Dialog(enable = not all_rear,tab="Inlet/outlet", group="Rear"));
  parameter Boolean all_floor = false
    "This flag enables all floor fluid connectors"
  annotation (Evaluate = true,Dialog(tab="Inlet/outlet", group="Floor"),choices(__Dymola_checkBox=true));
  parameter Integer ij_floor[:,2] = [0,0]
    "coordinates of the (floor) fluid connectors [i1,j1;i2,j2;...]"
  annotation (Evaluate = true,Dialog(enable = not all_floor,tab="Inlet/outlet", group="Floor"));
  parameter Boolean all_cei = false
    "This flag enables all ceiling fluid connectors"
  annotation (Evaluate = true,Dialog(tab="Inlet/outlet", group="Ceiling"),choices(__Dymola_checkBox=true));
  parameter Integer ij_cei[:,2] = [0,0]
    "coordinates of the (ceiling) fluid connectors [i1,j1;i2,j2;...]"
  annotation (Evaluate = true,Dialog(enable = not all_cei,tab="Inlet/outlet", group="Ceiling"));
  /*********** GEOMETRY ****************************************************************************/
  parameter Modelica.SIunits.Distance roomWidth(min=0)=3 "width of the room"
                                                               annotation (Dialog(tab="Geometry", group="sizes"));
  parameter Modelica.SIunits.Distance roomHeight(min=0)=2.5
    " height of the room"                                           annotation (Dialog(tab="Geometry", group="sizes"));
  parameter Modelica.SIunits.Distance roomBase(min=0)=3 "base of the room"
                                                              annotation (Dialog(tab="Geometry", group="sizes"));
  parameter Integer I(min=3) = 3 "number of volumes along x-direction" annotation (Dialog(tab="Geometry", group="grid"));
  parameter Integer J(min=3) = 3 "Number of volumes along y-direction" annotation (Dialog(tab="Geometry", group="grid"));
  parameter Integer K(min=3) = 3 "number of volumes along z-direction" annotation (Dialog(tab="Geometry", group="grid"));
  parameter Real X_frac[I] = 1/I*ones(I)
    "fraction of base assigned to each volume (the sum must be one)" annotation (Dialog(tab="Geometry", group="grid"));
  parameter Real Y_frac[J] = 1/J*ones(J)
    "fraction of width assigned to each volume (the sum must be one)" annotation (Dialog(tab="Geometry", group="grid"));
  parameter Real Z_frac[K] = 1/K*ones(K)
    "fraction of height assigned to each volume (the sum must be one)" annotation (Dialog(tab="Geometry", group="grid"));
  parameter Modelica.SIunits.Length h_ref=0
    "height of the bottom of the room with respect to a given reference" annotation (Dialog(tab="Geometry", group="reference height"));
  /************ INITIAL CONDITIONS +++**************************************************************/
  parameter Modelica.SIunits.Temperature Tstart(displayUnit="K")=273.15 + 20
    "Fluid initial temperature"  annotation (Dialog(tab="Initialisation", group="Fluid"));
  parameter Modelica.SIunits.Velocity Vxstart=0
    "initial fluid velocity (x-direction)"                              annotation (Dialog(tab="Initialisation", group="Fluid"));
  parameter Modelica.SIunits.Velocity Vystart=0
    "initial fluid velocity (y-direction)"                              annotation (Dialog(tab="Initialisation", group="Fluid"));
  parameter Modelica.SIunits.Velocity Vzstart=0
    "initial fluid velocity (z-direction)"                              annotation (Dialog(tab="Initialisation", group="Fluid"));
  /*************** TURBULENCE MODEL +***************************************************************/
  parameter Boolean laminar = false "flag that avoid the turbulence model"
                                           annotation (Dialog(tab="Turbulence", group="model"),choices(__Dymola_checkBox=true));
  parameter Real Kturb(min=0) = 0.03874 "Coefficient of th eturbulence model"
                                           annotation (Dialog(enable = not laminar, tab="Turbulence", group="model"));
  parameter String lengthWall = "Xu" "Method for the wall distance"
                                    annotation (Dialog(enable = not laminar, tab="Turbulence", group="model"),
                                    choices(choice=1 "Xu",choice=2 "Xu"));
  /***************** CD EQUATIONS ******************************************************************/
  parameter Boolean convective = true
    "Flag that enable the convective term in the CD eqns"
    annotation (Evaluate = true,Dialog(tab="CD equations", group="enabling terms"),choices(__Dymola_checkBox=true));
  parameter Boolean energy = true "Include the energy equation"
    annotation (Evaluate = true,Dialog(tab="CD equations", group="enabling terms"),choices(__Dymola_checkBox=true));
  parameter Boolean gravity = true "Include the gravity terms"
    annotation (Evaluate = true,Dialog(tab="CD equations", group="enabling terms"),choices(__Dymola_checkBox=true));
  parameter CFDModelica.Units.IntScheme intScheme=CFDModelica.Units.IntScheme.Hybrid
    "Method employed for integrating the convective term" annotation (Evaluate=
        true, Dialog(
      enable=convective,
      tab="CD equations",
      group="discretisation"));
  parameter CFDModelica.Units.IntScheme intSchemeEnergy=CFDModelica.Units.IntScheme.Hybrid
    "Method employed for integrating the convective term in the energy equation"
    annotation (Evaluate=true, Dialog(
      enable=convective and energy,
      tab="CD equations",
      group="discretisation"));
  /***************** BOUNDARIES *********************************************************************/
  parameter Real Kll(min=0) = 0.41
    "Coefficient for the Log Layer wall function"
     annotation (Evaluate = true,Dialog(enable = not laminar, tab="Boundaries", group="Wall function"));
  parameter Real E(min=3) = 8.9 "Coefficient for the Log Layer wall function"
     annotation (Evaluate = true,Dialog(enable = not laminar, tab="Boundaries", group="Wall function"));
  parameter Real y_plus_min(min=3) = 5 "Minimum y+"
     annotation (Evaluate = true,Dialog(enable = not laminar, tab="Boundaries", group="Wall function"));
  parameter Real Y_PLUS(min=3) = 30 "Treshold value for y+"
     annotation (Evaluate = true,Dialog(enable = not laminar, tab="Boundaries", group="Wall function"));
  parameter Boolean Hmean = false "Averaged value of the h wall"
     annotation (Evaluate = true,Dialog(enable = not Hfixed, tab="Boundaries", group="Averaged"),choices(__Dymola_checkBox=true));
  parameter Real T_h_mean(min=3) = 5 "Time constant"
     annotation (Evaluate = true,Dialog(enable = Hmean, tab="Boundaries", group="Averaged"));
  parameter Boolean Hfixed = true "Averaged value of the h wall"
     annotation (Evaluate = true,Dialog(enable = not Hmean, tab="Boundaries", group="Fixed"),choices(__Dymola_checkBox=true));
  parameter Real Hright(min=0) = 5
    "Fixed convective heat transfer coefficient (right)"
     annotation (Evaluate = true,Dialog(enable = (Hfixed), tab="Boundaries", group="Fixed"));
  parameter Real Hleft(min=0) = 5
    "Fixed convective heat transfer coefficient (left)"
     annotation (Evaluate = true,Dialog(enable = (Hfixed), tab="Boundaries", group="Fixed"));
  parameter Real Hfloor(min=0) = 5
    "Fixed convective heat transfer coefficient (floor)"
     annotation (Evaluate = true,Dialog(enable = (Hfixed), tab="Boundaries", group="Fixed"));
  parameter Real Hcei(min=0) = 5
    "Fixed convective heat transfer coefficient (ceiling)"
     annotation (Evaluate = true,Dialog(enable = (Hfixed), tab="Boundaries", group="Fixed"));
  parameter Real Hfront(min=0) = 5
    "Fixed convective heat transfer coefficient (front)"
     annotation (Evaluate = true,Dialog(enable = (Hfixed), tab="Boundaries", group="Fixed"));
  parameter Real Hrear(min=0) = 5
    "Fixed convective heat transfer coefficient (rear)"
     annotation (Evaluate = true,Dialog(enable = (Hfixed), tab="Boundaries", group="Fixed"));
  /************ OUTPUT SELECTION **********************************************************************/
  parameter Boolean useOUTPUTS = false
    "|Output preferences| This flag enables the output for post processing" annotation(Evaluate = true);
  Modelica.Blocks.Interfaces.RealOutput T_[I,J,K] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput rho_[I,J,K] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput Vx_[I+1,J+2,K+2] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput Vy_[I+2,J+1,K+2] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput Vz_[I+2,J+2,K+1] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput P_[I,J,K] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput I_ if  useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput J_ if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput K_ if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput roomBase_ if  useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput roomWidth_ if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput roomHeight_ if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput X_frac_[I] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput Y_frac_[J] if useOUTPUTS;
  Modelica.Blocks.Interfaces.RealOutput Z_frac_[K] if useOUTPUTS;
  /************ INTERFACES ****************************************************************************/
  CFDModelica.Interfaces.HeatPort heatPort_left[J,K]
    "matrix of heat port connectors, left wall"
    annotation (Placement(transformation(extent={{-100,10},{-80,30}})));
  CFDModelica.Interfaces.HeatPort heatPort_right[J,K]
    "matrix of heat port connectors, right wall"
    annotation (Placement(transformation(extent={{80,10},{100,30}})));
  CFDModelica.Interfaces.HeatPort heatPort_floor[I,J]
    "matrix of heat port connectors, floor"
    annotation (Placement(transformation(extent={{-30,-100},{-10,-80}})));
  CFDModelica.Interfaces.HeatPort heatPort_cei[I,J]
    "matrix of heat port connectors, ceiling"
    annotation (Placement(transformation(extent={{-30,80},{-10,100}})));
  CFDModelica.Interfaces.HeatPort heatPort_rear[I,K]
    "matrix of heat port connectors, rear wall"
    annotation (Placement(transformation(extent={{60,80},{80,100}})));
  CFDModelica.Interfaces.HeatPort heatPort_front[I,K]
    "matrix of heat port connectors, front wall"
    annotation (Placement(transformation(extent={{-100,-80},{-80,-60}})));
  CFDModelica.Interfaces.HeatPort heatPort[I,J,K]
    "matrix of heat port connectors, connected to each volume (heat sources)"
    annotation (Placement(transformation(extent={{80,-100},{100,-80}})));
  CFDModelica.Interfaces.flangeA fluidPort_left[J,K](redeclare package Medium
      = Medium) "matrix of fluid connectors, left wall"
    annotation (Placement(transformation(extent={{-100,-30},{-80,-10}})));
  CFDModelica.Interfaces.flangeB fluidPort_right[J,K](redeclare package Medium
      = Medium) "matrix of fluid connectors, right wall"
    annotation (Placement(transformation(extent={{80,-30},{100,-10}})));
  CFDModelica.Interfaces.flangeA fluidPort_floor[I,J](redeclare package Medium
      = Medium) "matrix of fluid connectors, floor wall"
    annotation (Placement(transformation(extent={{10,-100},{30,-80}})));
  CFDModelica.Interfaces.flangeB fluidPort_cei[I,J](redeclare package Medium =
        Medium) "matrix of fluid connectors, ceiling wall"
    annotation (Placement(transformation(extent={{10,80},{30,100}})));
  CFDModelica.Interfaces.flangeA fluidPort_front[I,K](redeclare package Medium
      = Medium) "matrix of fluid connectors, front wall"
    annotation (Placement(transformation(extent={{-80,-100},{-60,-80}})));
  CFDModelica.Interfaces.flangeB fluidPort_rear[I,K](redeclare package Medium
      = Medium) "matrix of fluid connectors, rear wall"
    annotation (Placement(transformation(extent={{80,60},{100,80}})));
  /***************** VARIABLES **********************************************************************/
  Modelica.SIunits.Mass Mtot "mass of air contained within the room";
  Modelica.SIunits.Energy Etot "Energy contained within the room";
  Modelica.SIunits.Power FluxLeft "heat flow rate fron the left wall";
  Modelica.SIunits.Power FluxRight "heat flow rate fron the right wall";
  Modelica.SIunits.Power FluxFloor "heat flow rate fron the floor";
  Modelica.SIunits.Power FluxCeiling "heat flow rate fron the ceiling";
  Modelica.SIunits.Power FluxFront "heat flow rate fron the front wall";
  Modelica.SIunits.Power FluxRear "heat flow rate fron the rear wall";
  /***************** FLUID PROPERTIES *********************************************************************/
  Modelica.SIunits.Temperature T[I + 2,J + 2,K + 2](each start=Tstart, each
      stateSelect=StateSelect.prefer)
    "temperatures inside the volume and boundaries"                                                  annotation(displayUnit="K");
  Modelica.SIunits.Density rho[I + 2,J + 2,K + 2]
    "densities inside the volume and boundaries ";
  Modelica.SIunits.Pressure P[I,J,K](displayUnit="Pa", each stateSelect=
        StateSelect.prefer) "pressures inside the volume";
  Modelica.SIunits.Mass m[I,J,K] "masses within the volume";
  Modelica.SIunits.SpecificEnergy e[I,J,K] "specific energy inside the volume";
  Modelica.SIunits.Velocity Vx[I + 1,J + 2,K + 2](each start=Vxstart)
    "fluid velocities inside the volume and boundaries";
  Modelica.SIunits.Velocity Vy[I + 2,J + 1,K + 2](each start=Vystart)
    "fluid velocities inside the volume and boundaries";
  Modelica.SIunits.Velocity Vz[I + 2,J + 2,K + 1](each start=Vzstart)
    "fluid velocities inside the volume and boundaries";
  Modelica.SIunits.DynamicViscosity MUx[I + 1,J + 2,K + 2](each start=mu)
    "fluid turbulence viscosity inside the volume and boundaries";
  Modelica.SIunits.DynamicViscosity MUy[I + 2,J + 1,K + 2](each start=mu)
    "fluid turbulence viscosity inside the volume and boundaries";
  Modelica.SIunits.DynamicViscosity MUz[I + 2,J + 2,K + 1](each start=mu)
    "fluid turbulence viscosity inside the volume and boundaries";
  /***************** CONVECTIVE HEAT TRANSFER COEFFICIENTS ***********************************************/
  Modelica.SIunits.CoefficientOfHeatTransfer h_left[J,K]
    "matrix of convetive heat transfer coefficients";
  Modelica.SIunits.CoefficientOfHeatTransfer h_right[J,K]
    "matrix of convetive heat transfer coefficients";
  Modelica.SIunits.CoefficientOfHeatTransfer h_cei[I,J]
    "matrix of convetive heat transfer coefficients";
  Modelica.SIunits.CoefficientOfHeatTransfer h_floor[I,J]
    "matrix of convetive heat transfer coefficients";
  Modelica.SIunits.CoefficientOfHeatTransfer h_front[I,K]
    "matrix of convetive heat transfer coefficients";
  Modelica.SIunits.CoefficientOfHeatTransfer h_rear[I,K]
    "matrix of convetive heat transfer coefficients";
  Real h_avg_left "averaged convetive heat transfer coefficient";
  Real h_avg_right "averaged convetive heat transfer coefficient";
  Real h_avg_floor "averaged convetive heat transfer coefficient";
  Real h_avg_cei "averaged convetive heat transfer coefficient";
  Real h_avg_front "averaged convetive heat transfer coefficient";
  Real h_avg_rear "averaged convetive heat transfer coefficient";
protected
  Modelica.Blocks.Interfaces.RealOutput T_a[I,J,K];
  Modelica.Blocks.Interfaces.RealOutput rho_a[I,J,K];
  Modelica.Blocks.Interfaces.RealOutput Vx_a[I+1,J+2,K+2];
  Modelica.Blocks.Interfaces.RealOutput Vy_a[I+2,J+1,K+2];
  Modelica.Blocks.Interfaces.RealOutput Vz_a[I+2,J+2,K+1];
  Modelica.Blocks.Interfaces.RealOutput P_a[I,J,K];
  Modelica.Blocks.Interfaces.RealOutput I_a;
  Modelica.Blocks.Interfaces.RealOutput J_a;
  Modelica.Blocks.Interfaces.RealOutput K_a;
  Modelica.Blocks.Interfaces.RealOutput roomBase_a;
  Modelica.Blocks.Interfaces.RealOutput roomWidth_a;
  Modelica.Blocks.Interfaces.RealOutput roomHeight_a;
  Modelica.Blocks.Interfaces.RealOutput X_frac_a[I];
  Modelica.Blocks.Interfaces.RealOutput Y_frac_a[J];
  Modelica.Blocks.Interfaces.RealOutput Z_frac_a[K];

  /************** VARIABLES for connection sets/boundaries *************************************/
  // Real Tmix_left[J,K];
  // Real Tmix_right[J,K];
  // Real Tmix_front[I,K];
  // Real Tmix_rear[I,K];
  // Real Tmix_cei[I,J];
  // Real Tmix_floor[I,J];
  parameter CFDModelica.Units.NormPower alpha=1 annotation (Evaluate=true);
  parameter CFDModelica.Units.NormMflowRate beta=1 annotation (Evaluate=true);
  Modelica.SIunits.Pressure Pleft[J,K](displayUnit="Pa");
  Modelica.SIunits.Pressure Pright[J,K](displayUnit="Pa");
  Modelica.SIunits.Pressure Pfront[I,K](displayUnit="Pa");
  Modelica.SIunits.Pressure Prear[I,K](displayUnit="Pa");
  Modelica.SIunits.Pressure Pfloor[I,J](displayUnit="Pa");
  Modelica.SIunits.Pressure Pcei[I,J](displayUnit="Pa");
  Real De_bl_x[J,K];
  Real Fe_bl_x[J,K];
  Real Pe_bl_x[J,K];
  Real aE_bl_x[J,K];
  Real Dw_br_x[J,K];
  Real Fw_br_x[J,K];
  Real Pw_br_x[J,K];
  Real aW_br_x[J,K];

  Real Dn_bf_y[I,K];
  Real Fn_bf_y[I,K];
  Real Pn_bf_y[I,K];
  Real aN_bf_y[I,K];
  Real Ds_bre_y[I,K];
  Real Fs_bre_y[I,K];
  Real Ps_bre_y[I,K];
  Real aS_bre_y[I,K];

  Real Dt_bb_z[I,J];
  Real Ft_bb_z[I,J];
  Real Pt_bb_z[I,J];
  Real aT_bb_z[I,J];
  Real Db_bt_z[I,J];
  Real Fb_bt_z[I,J];
  Real Pb_bt_z[I,J];
  Real aB_bt_z[I,J];

  /***************** VARIABLES FOR MASS PRESERVATION EQUATIONS *********************************************************************/
  /***************** Diffusive conductances, Peclet Numbers, Mass Flows and A-coefficients *****************************************/

  Real De_M[I,J,K];
  Real Dw_M[I,J,K];
  Real Ds_M[I,J,K];
  Real Dn_M[I,J,K];
  Real Db_M[I,J,K];
  Real Dt_M[I,J,K];

  Real Pe_M[I,J,K];
  Real Pw_M[I,J,K];
  Real Ps_M[I,J,K];
  Real Pn_M[I,J,K];
  Real Pb_M[I,J,K];
  Real Pt_M[I,J,K];

  Real Fe_M[I,J,K];
  Real Fw_M[I,J,K];
  Real Fs_M[I,J,K];
  Real Fn_M[I,J,K];
  Real Fb_M[I,J,K];
  Real Ft_M[I,J,K];

  Real aE_M[I,J,K];
  Real aW_M[I,J,K];
  Real aB_M[I,J,K];
  Real aT_M[I,J,K];
  Real aS_M[I,J,K];
  Real aN_M[I,J,K];
  Real aP_M[I,J,K];
  /***************** VARIABLES FOR x-MOMENTUM EQUATIONS ****************************************************************************/
  /***************** Diffusive conductances, Peclet Numbers, Mass Flows and A-coefficients *****************************************/
  Real De_x[I-1,J,K];
  Real Dw_x[I-1,J,K];
  Real Ds_x[I-1,J,K];
  Real Dn_x[I-1,J,K];
  Real Db_x[I-1,J,K];
  Real Dt_x[I-1,J,K];
  Real Pe_x[I-1,J,K];
  Real Pw_x[I-1,J,K];
  Real Ps_x[I-1,J,K];
  Real Pn_x[I-1,J,K];
  Real Pb_x[I-1,J,K];
  Real Pt_x[I-1,J,K];

  Real Fe_x[I-1,J,K];
  Real Fw_x[I-1,J,K];
  Real Fs_x[I-1,J,K];
  Real Fn_x[I-1,J,K];
  Real Fb_x[I-1,J,K];
  Real Ft_x[I-1,J,K];

  Real aE_x[I-1,J,K];
  Real aW_x[I-1,J,K];
  Real aS_x[I-1,J,K];
  Real aN_x[I-1,J,K];
  Real aB_x[I-1,J,K];
  Real aT_x[I-1,J,K];
  Real aP_x[I-1,J,K];
  /***************** VARIABLES FOR y-MOMENTUM EQUATIONS ****************************************************************************/
  /***************** Diffusive conductances, Peclet Numbers, Mass Flows and A-coefficients *****************************************/
  Real De_y[I,J-1,K];
  Real Dw_y[I,J-1,K];
  Real Ds_y[I,J-1,K];
  Real Dn_y[I,J-1,K];
  Real Db_y[I,J-1,K];
  Real Dt_y[I,J-1,K];
  Real Pe_y[I,J-1,K];
  Real Pw_y[I,J-1,K];
  Real Ps_y[I,J-1,K];
  Real Pn_y[I,J-1,K];
  Real Pb_y[I,J-1,K];
  Real Pt_y[I,J-1,K];

  Real Fe_y[I,J-1,K];
  Real Fw_y[I,J-1,K];
  Real Fs_y[I,J-1,K];
  Real Fn_y[I,J-1,K];
  Real Fb_y[I,J-1,K];
  Real Ft_y[I,J-1,K];

  Real aE_y[I,J-1,K];
  Real aW_y[I,J-1,K];
  Real aS_y[I,J-1,K];
  Real aN_y[I,J-1,K];
  Real aB_y[I,J-1,K];
  Real aT_y[I,J-1,K];
  Real aP_y[I,J-1,K];
  /***************** VARIABLES FOR z-MOMENTUM EQUATIONS ****************************************************************************/
  /***************** Diffusive conductances, Peclet Numbers, Mass Flows and A-coefficients *****************************************/
  Real De_z[I,J,K-1];
  Real Dw_z[I,J,K-1];
  Real Ds_z[I,J,K-1];
  Real Dn_z[I,J,K-1];
  Real Db_z[I,J,K-1];
  Real Dt_z[I,J,K-1];
  Real Pe_z[I,J,K-1];
  Real Pw_z[I,J,K-1];
  Real Ps_z[I,J,K-1];
  Real Pn_z[I,J,K-1];
  Real Pb_z[I,J,K-1];
  Real Pt_z[I,J,K-1];

  Real Fe_z[I,J,K-1];
  Real Fw_z[I,J,K-1];
  Real Fs_z[I,J,K-1];
  Real Fn_z[I,J,K-1];
  Real Fb_z[I,J,K-1];
  Real Ft_z[I,J,K-1];

  Real aE_z[I,J,K-1];
  Real aW_z[I,J,K-1];
  Real aS_z[I,J,K-1];
  Real aN_z[I,J,K-1];
  Real aB_z[I,J,K-1];
  Real aT_z[I,J,K-1];
  Real aP_z[I,J,K-1];
  /***************** FLUID PROPERTIES ************************************************************************/
  parameter Modelica.SIunits.DynamicViscosity mu=Medium.eta_const
    "laminar fluid viscosity"                                                annotation(Evaluate = true);
  parameter Modelica.SIunits.ThermalConductivity gamma=Medium.lambda_const
    "fluid thermal conductivity" annotation(Evaluate = true);
  parameter Modelica.SIunits.Density rho_o=Medium.d0
    "fluid density linearisation value"                                   annotation(Evaluate = true);
  parameter Modelica.SIunits.SpecificHeatCapacity cv=Medium.cv_const
    "fluid specific heat at constant volume" annotation(Evaluate = true);
  parameter Real Pr = mu/gamma "prandtl number" annotation(Evaluate = true);
  /***************** GEOMETRIC PROPERTIES *******************************************************************/
  parameter Modelica.SIunits.Length height[I,J,K]=
      CFDModelica.Functions.CubeHeights(
      I,
      J,
      K,
      dz[:],
      h_ref) "height above the floor of each volume"
                                            annotation(Evaluate = true);
  parameter Modelica.SIunits.Acceleration g=if gravity then Modelica.Constants.g_n
       else 0 "constant gravity acceleration"
                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance dx[I]=roomBase*X_frac[:];
  parameter Modelica.SIunits.Distance dy[J]=roomWidth*Y_frac[:];
  parameter Modelica.SIunits.Distance dz[K]=roomHeight*Z_frac[:];
  parameter Modelica.SIunits.Distance dzx[K + 1]=CFDModelica.Functions.Dzx( K,
      dz[:]) "distance between x-velocity along z direction"
                                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance dzy[K + 1]=CFDModelica.Functions.Dzy( K,
      dz[:]) "distance between y-velocity along z direction"
                                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance dyx[J + 1]=CFDModelica.Functions.Dyx( J,
      dy[:]) "distance between x-velocity along y direction"
                                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance dxy[I + 1]=CFDModelica.Functions.Dxy( I,
      dx[:]) "distance between y-velocity along x direction"
                                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance dxz[I + 1]=CFDModelica.Functions.Dxz( I,
      dx[:]) "distance between z-velocity along x direction"
                                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance dyz[J + 1]=CFDModelica.Functions.Dyz( J,
      dy[:]) "distance between z-velocity along y direction"
                                                    annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance y_left=dx[1]/2
    "distance for the boundary layer (left)";
  parameter Modelica.SIunits.Distance y_right=dx[I]/2
    "distance for the boundary layer (right)" annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance y_front=dy[1]/2
    "distance for the boundary layer (front)" annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance y_rear=dy[J]/2
    "distance for the boundary layer (rear)";
  parameter Modelica.SIunits.Distance y_floor=dz[1]/2
    "distance for the boundary layer (floor)" annotation(Evaluate = true);
  parameter Modelica.SIunits.Distance y_cei=dz[K]/2
    "distance for the boundary layer (ceiling)" annotation(Evaluate = true);
  annotation (Diagram(graphics), Icon(graphics));
end Room3D_CFDwall_base;
