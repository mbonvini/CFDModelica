within CFDModelica.Flow3D;
partial model Room3D_CFDwall_massANDenergy
  //extends OOMS_BEPI.SubZonalModels_Extension.Flow3D_a.Room3D_CFDwall_walls;
  extends Room3D_CFDwall_base;
equation
  for i in 2:I+1 loop
    for j in 2:J+1 loop
      for k in 2:K + 1 loop
        // mass flows
        Fe_M[i - 1, j - 1, k - 1] = dy[j - 1]*dz[k
           - 1]*Vx[i, j, k]*(if noEvent(Vx[i, j, k]
           > 0) then rho[i, j, k] else rho[i + 1, j,
          k]);
        Fw_M[i - 1, j - 1, k - 1] = dy[j - 1]*dz[k
           - 1]*Vx[i - 1, j, k]*(if noEvent(Vx[i - 1, j,
          k] > 0) then rho[i - 1, j, k]
           else rho[i, j, k]);
        Fn_M[i - 1, j - 1, k - 1] = dz[k -
          1]*dx[i - 1]*Vy[i, j, k]*(if noEvent(Vy[i, j,
          k] > 0) then rho[i, j, k] else
          rho[i, j + 1, k]);
        Fs_M[i - 1, j - 1, k - 1] = dz[k -
          1]*dx[i - 1]*Vy[i, j - 1, k]*(if noEvent(Vy[i, j -
          1, k] > 0) then rho[i, j - 1, k]
           else rho[i, j, k]);
        Ft_M[i - 1, j - 1, k - 1] = dy[j - 1]*dx[i - 1]*Vz[i,
          j, k]*(if noEvent(Vz[i, j, k] >
          0) then rho[i, j, k] else rho[i, j, k
           + 1]);
        Fb_M[i - 1, j - 1, k - 1] = dy[j - 1]*dx[i - 1]*Vz[i,
          j, k - 1]*(if noEvent(Vz[i, j, k
           - 1] > 0) then rho[i, j, k - 1] else rho[i, j,
          k]);

        // Peclet numbers
        Pe_M[i - 1, j - 1, k - 1] = if noEvent(abs(De_M[i - 1,
          j - 1, k - 1]) < Modelica.Constants.eps) then 1
           else Fe_M[i - 1, j - 1, k - 1]/De_M[i - 1, j - 1,
          k - 1];
        Pw_M[i - 1, j - 1, k - 1] = if noEvent(abs(Dw_M[i - 1,
          j - 1, k - 1]) < Modelica.Constants.eps) then 1
           else Fw_M[i - 1, j - 1, k - 1]/Dw_M[i - 1, j - 1,
          k - 1];
        Pn_M[i - 1, j - 1, k - 1] = if noEvent(abs(Dn_M[i - 1,
          j - 1, k - 1]) < Modelica.Constants.eps) then 1
           else Fn_M[i - 1, j - 1, k - 1]/Dn_M[i - 1, j - 1,
          k - 1];
        Ps_M[i - 1, j - 1, k - 1] = if noEvent(abs(Ds_M[i - 1,
          j - 1, k - 1]) < Modelica.Constants.eps) then 1
           else Fs_M[i - 1, j - 1, k - 1]/Ds_M[i - 1, j - 1,
          k - 1];
        Pt_M[i - 1, j - 1, k - 1] = if noEvent(abs(Dt_M[i - 1,
          j - 1, k - 1]) < Modelica.Constants.eps) then 1
           else Ft_M[i - 1, j - 1, k - 1]/Dt_M[i - 1, j - 1,
          k - 1];
        Pb_M[i - 1, j - 1, k - 1] = if noEvent(abs(Db_M[i - 1,
          j - 1, k - 1]) < Modelica.Constants.eps) then 1
           else Fb_M[i - 1, j - 1, k - 1]/Db_M[i - 1, j - 1,
          k - 1];

        // terms of the energy preservation equation
        aE_M[i - 1, j - 1, k - 1] = De_M[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A((Medium.specificHeatCapacityCp(
          Medium.ThermodynamicState(P[i - 1, j - 1, k - 1], T[i, j, k])))*Pe_M[
          i - 1, j - 1, k - 1], intSchemeEnergy) + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j -
          1, k - 1], T[i, j, k])))*max(-Fe_M[i - 1, j - 1, k - 1], 0);

        aW_M[i - 1, j - 1, k - 1] = Dw_M[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A((Medium.specificHeatCapacityCp(
          Medium.ThermodynamicState(P[i - 1, j - 1, k - 1], T[i, j, k])))*Pw_M[
          i - 1, j - 1, k - 1], intSchemeEnergy) + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j -
          1, k - 1], T[i, j, k])))*max(Fw_M[i - 1, j - 1, k - 1], 0);

        aN_M[i - 1, j - 1, k - 1] = Dn_M[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A((Medium.specificHeatCapacityCp(
          Medium.ThermodynamicState(P[i - 1, j - 1, k - 1], T[i, j, k])))*Pn_M[
          i - 1, j - 1, k - 1], intSchemeEnergy) + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j -
          1, k - 1], T[i, j, k])))*max(-Fn_M[i - 1, j - 1, k - 1], 0);

        aS_M[i - 1, j - 1, k - 1] = Ds_M[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A((Medium.specificHeatCapacityCp(
          Medium.ThermodynamicState(P[i - 1, j - 1, k - 1], T[i, j, k])))*Ps_M[
          i - 1, j - 1, k - 1], intSchemeEnergy) + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j -
          1, k - 1], T[i, j, k])))*max(Fs_M[i - 1, j - 1, k - 1], 0);

        aT_M[i - 1, j - 1, k - 1] = Dt_M[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A((Medium.specificHeatCapacityCp(
          Medium.ThermodynamicState(P[i - 1, j - 1, k - 1], T[i, j, k])))*Pt_M[
          i - 1, j - 1, k - 1], intSchemeEnergy) + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j -
          1, k - 1], T[i, j, k])))*max(-Ft_M[i - 1, j - 1, k - 1], 0);

        aB_M[i - 1, j - 1, k - 1] = Db_M[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A((Medium.specificHeatCapacityCp(
          Medium.ThermodynamicState(P[i - 1, j - 1, k - 1], T[i, j, k])))*Pb_M[
          i - 1, j - 1, k - 1], intSchemeEnergy) + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j -
          1, k - 1], T[i, j, k])))*max(Fb_M[i - 1, j - 1, k - 1], 0);

        aP_M[i - 1, j - 1, k - 1] = aE_M[i - 1, j - 1,
          k - 1] + aW_M[i - 1, j - 1, k -
          1] + aN_M[i - 1, j - 1, k - 1] + aS_M[i - 1, j - 1,
          k - 1] + aB_M[i - 1, j - 1, k -
          1] + aT_M[i - 1, j - 1, k - 1] + (
          Medium.specificHeatCapacityCp(Medium.ThermodynamicState(P[i - 1, j - 1,
          k - 1], T[i, j, k])))*(Fe_M[i -
          1, j - 1, k - 1] - Fw_M[i - 1, j - 1, k
           - 1] + Fn_M[i - 1, j - 1, k - 1] - Fs_M[i - 1, j -
          1, k - 1] + Ft_M[i - 1, j - 1, k
           - 1] - Fb_M[i - 1, j - 1, k - 1]);

        // The PV=nRT relationship (valid for ideal gases) had ben linearised by P,T, obtaining:
        //rho[i,j,k] = rho_o + ComprCoeff*P[i-1,j-1,k-1] - ThermalExpCoeff*(T[i,j,k]-To);
        //rho[i,j,k] = if energy then Medium.density_pT(P[i-1,j-1,k-1],T[i,j,k]) else Medium.density_pT(P[i-1,j-1,k-1],Medium.T0);
        rho[i, j, k] = Medium.density_pT(P[i - 1, j - 1, k - 1], T[i, j, k]);

        // mass of the air volume
        m[i - 1, j - 1, k - 1] = rho[i, j, k]
          *dx[i - 1]*dz[k - 1]*dy[j - 1];

        // mass conservation
        der(m[i - 1, j - 1, k - 1]) = Fw_M[i - 1, j - 1,
          k - 1] - Fe_M[i - 1, j - 1, k -
          1] + Fs_M[i - 1, j - 1, k - 1] - Fn_M[i - 1, j - 1,
          k - 1] + Fb_M[i - 1, j - 1, k -
          1] - Ft_M[i - 1, j - 1, k - 1];

        // specific energy
        //e[i - 1, j - 1, k - 1] = cv*T[i, j, k];
        e[i - 1, j - 1, k - 1] =
          Medium.specificInternalEnergy_pT(P[i - 1, j - 1, k -
          1], T[i, j, k]);

        // energy balance
        // 1) rigorous expression:
        //     der(m[i-1,J,k-1]*e[i-1,J,k-1])
        // 2) simplified expression:
        //     rho_o*cv*dx[i-1]*dz[k-1]*dy[j-1]*der(T[i,j,k])
        rho_o*dx[i - 1]*dz[k - 1]*dy[j - 1]*der(e[i - 1, j -
          1, k - 1]) = if energy then -aP_M[i - 1, j - 1,
          k - 1]*T[i, j, k] + aE_M[i - 1,
          j - 1, k - 1]*T[i + 1, j, k] +
          aW_M[i - 1, j - 1, k - 1]*T[i - 1, j, k]
           + aN_M[i - 1, j - 1, k - 1]*T[i, j + 1, k]
           + aS_M[i - 1, j - 1, k - 1]*T[i, j - 1, k]
           + aT_M[i - 1, j - 1, k - 1]*T[i, j, k
           + 1] + aB_M[i - 1, j - 1, k - 1]*T[i, j, k
           - 1] + heatPort[i - 1, j - 1, k - 1].Q_flow else 0;
         /*
         rho_o*cv*dx[i-1]*dy[j-1]*dz[k-1]*der(T[i,j,k]) =
            if energy then 
              - aP_M[i-1,j-1,k-1]*T[i,j,k]
              + aE_M[i-1,j-1,k-1]*T[i+1,j,k]
              + aW_M[i-1,j-1,k-1]*T[i-1,j,k]
              + aN_M[i-1,j-1,k-1]*T[i,j+1,k]
              + aS_M[i-1,j-1,k-1]*T[i,j-1,k]
              + aT_M[i-1,j-1,k-1]*T[i,j,k+1]
              + aB_M[i-1,j-1,k-1]*T[i,j,k-1]
              + heatPort[i-1,j-1,k-1].Q_flow else 
              0;
         */

         // thermal connectors
        heatPort[i - 1, j - 1, k - 1].T = T[i, j, k];

      end for;
    end for;
  end for;

end Room3D_CFDwall_massANDenergy;
