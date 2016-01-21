within CFDModelica.Flow3D;
partial model Room3D_CFDwall_momentum
  // extends OOMS_BEPI.SubZonalModels_Extension.Flow3D_a.Room3D_CFDwall_turbulence;
  extends Room3D_CFDwall_base;
equation
  ////////////////////////////////////////
  // X-MOMENTUM PRESERVATION
  for i in 2:I loop
    for j in 2:J+1 loop
      for k in 2:K + 1 loop
        // diffusive fluxes
        De_x[i - 1, j - 1, k - 1] = dy[j - 1]*dz[k
           - 1]*(0.5*MUx[i, j, k] + 0.5*MUx[i + 1, j,
          k])/dx[i];
        Dw_x[i - 1, j - 1, k - 1] = dy[j - 1]*dz[k
           - 1]*(0.5*MUx[i, j, k] + 0.5*MUx[i - 1, j,
          k])/dx[i - 1];
        Dn_x[i - 1, j - 1, k - 1] = dz[k -
          1]*(0.5*dx[i] + 0.5*dx[i - 1])*(0.5*MUx[i, j, k] +
          0.5*MUx[i, j + 1, k])/dyx[j];
        Ds_x[i - 1, j - 1, k - 1] = dz[k -
          1]*(0.5*dx[i] + 0.5*dx[i - 1])*(0.5*MUx[i, j, k] +
          0.5*MUx[i, j - 1, k])/dyx[j - 1];
        Dt_x[i - 1, j - 1, k - 1] = dy[j - 1]*(0.5*dx[i] + 0.5
          *dx[i - 1])*(0.5*MUx[i, j, k] + 0.5*MUx[i, j,
          k + 1])/dzx[k];
        Db_x[i - 1, j - 1, k - 1] = dy[j - 1]*(0.5*dx[i] + 0.5
          *dx[i - 1])*(0.5*MUx[i, j, k] + 0.5*MUx[i, j,
          k - 1])/dzx[k - 1];

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fe_x[i - 1, j - 1, k - 1] = if convective then (Fw_M[i, j - 1,
          k - 1] + Fe_M[i, j - 1, k - 1])/
          2 else 0;
        Fw_x[i - 1, j - 1, k - 1] = if convective then (Fw_M[
          i - 1, j - 1, k - 1] + Fe_M[i - 1, j - 1, k
           - 1])/2 else 0;
        Fn_x[i - 1, j - 1, k - 1] = if convective then (dx[i -
          1]*Fn_M[i - 1, j - 1, k - 1] + dx[i]*Fn_M[i, j - 1,
          k - 1])/(dx[i - 1] + dx[i]) else 0;
        Fs_x[i - 1, j - 1, k - 1] = if convective then (dx[i -
          1]*Fs_M[i - 1, j - 1, k - 1] + dx[i]*Fs_M[i, j - 1,
          k - 1])/(dx[i - 1] + dx[i]) else 0;
        Ft_x[i - 1, j - 1, k - 1] = if convective then (dx[i -
          1]*Ft_M[i - 1, j - 1, k - 1] + dx[i]*Ft_M[i, j - 1,
          k - 1])/(dx[i - 1] + dx[i]) else 0;
        Fb_x[i - 1, j - 1, k - 1] = if convective then (dx[i -
          1]*Fb_M[i - 1, j - 1, k - 1] + dx[i]*Fb_M[i, j - 1,
          k - 1])/(dx[i - 1] + dx[i]) else 0;

        // Peclet number
        Pe_x[i - 1, j - 1, k - 1] = Fe_x[i - 1, j - 1,
          k - 1]/De_x[i - 1, j - 1, k - 1];
        Pw_x[i - 1, j - 1, k - 1] = Fw_x[i - 1, j - 1,
          k - 1]/Dw_x[i - 1, j - 1, k - 1];
        Pn_x[i - 1, j - 1, k - 1] = Fn_x[i - 1, j - 1,
          k - 1]/Dn_x[i - 1, j - 1, k - 1];
        Ps_x[i - 1, j - 1, k - 1] = Fs_x[i - 1, j - 1,
          k - 1]/Ds_x[i - 1, j - 1, k - 1];
        Pt_x[i - 1, j - 1, k - 1] = Ft_x[i - 1, j - 1,
          k - 1]/Dt_x[i - 1, j - 1, k - 1];
        Pb_x[i - 1, j - 1, k - 1] = Fb_x[i - 1, j - 1,
          k - 1]/Db_x[i - 1, j - 1, k - 1];

        // terms of the momentum equation
        aE_x[i - 1, j - 1, k - 1] = De_x[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pe_x[i - 1, j - 1, k - 1], intScheme) + max(-
          Fe_x[i - 1, j - 1, k - 1], 0);
        aW_x[i - 1, j - 1, k - 1] = Dw_x[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pw_x[i - 1, j - 1, k - 1], intScheme) + max(
          Fw_x[i - 1, j - 1, k - 1], 0);
        aN_x[i - 1, j - 1, k - 1] = Dn_x[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pn_x[i - 1, j - 1, k - 1], intScheme) + max(-
          Fn_x[i - 1, j - 1, k - 1], 0);
        aS_x[i - 1, j - 1, k - 1] = Ds_x[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Ps_x[i - 1, j - 1, k - 1], intScheme) + max(
          Fs_x[i - 1, j - 1, k - 1], 0);
        aT_x[i - 1, j - 1, k - 1] = Dt_x[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pt_x[i - 1, j - 1, k - 1], intScheme) + max(-
          Ft_x[i - 1, j - 1, k - 1], 0);
        aB_x[i - 1, j - 1, k - 1] = Db_x[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pb_x[i - 1, j - 1, k - 1], intScheme) + max(
          Fb_x[i - 1, j - 1, k - 1], 0);

        // Here the terms
        // + Fe_x[i,J,k] - Fw_x[i,J,k] + Fn_x[i,J,k] - Fs_x[i,J,k] + Ft_x[i,J,k] - Fb_x[i,J,k]
        // to be added are disregarded since the mass variation is not taken into account
        aP_x[i - 1, j - 1, k - 1] = aE_x[i - 1, j - 1, k - 1]
           + aW_x[i - 1, j - 1, k - 1] + aN_x[i - 1, j - 1,
          k - 1] + aS_x[i - 1, j - 1, k -
          1] + aB_x[i - 1, j - 1, k - 1] + aT_x[i - 1, j - 1,
          k - 1];

        // x-momentum equation
        (dx[i - 1]*0.5 + dx[i]*0.5)*dy[j - 1]*dz[k - 1]*
          rho_o*der(Vx[i, j, k]) + aP_x[i - 1, j - 1,
          k - 1]*Vx[i, j, k] = aE_x[i - 1,
          j - 1, k - 1]*Vx[i + 1, j, k] +
          aW_x[i - 1, j - 1, k - 1]*Vx[i - 1, j, k]
           + aN_x[i - 1, j - 1, k - 1]*Vx[i, j + 1, k]
           + aS_x[i - 1, j - 1, k - 1]*Vx[i, j - 1, k]
           + aT_x[i - 1, j - 1, k - 1]*Vx[i, j, k
           + 1] + aB_x[i - 1, j - 1, k - 1]*Vx[i, j,
          k - 1] + (P[i - 1, j - 1, k - 1]
           - P[i, j - 1, k - 1])*dy[j - 1]*dz[k
           - 1];

      end for;
    end for;
  end for;
  // END X-MOMENTUM PRESERVATION
  ////////////////////////////////////////

  ////////////////////////////////////////
  // Y-MOMENTUM PRESERVATION
  for i in 2:I+1 loop
    for j in 2:J loop
      for k in 2:K + 1 loop
        // diffusive fluxes
        De_y[i - 1, j - 1, k - 1] = dz[k -
          1]*(dy[j - 1]*0.5 + dy[j]*0.5)*(0.5*MUy[i, j, k] +
          0.5*MUy[i + 1, j, k])/dxy[i];
        Dw_y[i - 1, j - 1, k - 1] = dz[k -
          1]*(dy[j - 1]*0.5 + dy[j]*0.5)*(0.5*MUy[i, j, k] +
          0.5*MUy[i - 1, j, k])/dxy[i - 1];
        Dn_y[i - 1, j - 1, k - 1] = dx[i - 1]*dz[k
           - 1]*(0.5*MUy[i, j, k] + 0.5*MUy[i, j + 1,
          k])/dy[j];
        Ds_y[i - 1, j - 1, k - 1] = dx[i - 1]*dz[k
           - 1]*(0.5*MUy[i, j, k] + 0.5*MUy[i, j - 1,
          k])/dy[j - 1];
        Dt_y[i - 1, j - 1, k - 1] = dx[i - 1]*(dy[j - 1]*0.5 +
          dy[j]*0.5)*(0.5*MUy[i, j, k] + 0.5*MUy[i, j,
          k + 1])/dzy[k];
        Db_y[i - 1, j - 1, k - 1] = dx[i - 1]*(dy[j - 1]*0.5 +
          dy[j]*0.5)*(0.5*MUy[i, j, k] + 0.5*MUy[i, j,
          k - 1])/dzy[k - 1];

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fe_y[i - 1, j - 1, k - 1] = if convective then (dy[j - 1]*Fe_M[i
           - 1, j - 1, k - 1] + dy[j]*Fe_M[i - 1, j,
          k - 1])/(dy[j - 1] + dy[j]) else 0;
        Fw_y[i - 1, j - 1, k - 1] = if convective then (dy[j -
          1]*Fw_M[i - 1, j - 1, k - 1] + dy[j]*Fw_M[i - 1, j,
          k - 1])/(dy[j - 1] + dy[j]) else 0;
        Fn_y[i - 1, j - 1, k - 1] = if convective then (Fn_M[
          i - 1, j, k - 1] + Fs_M[i - 1, j, k
           - 1])/2 else 0;
        Fs_y[i - 1, j - 1, k - 1] = if convective then (Fn_M[
          i - 1, j - 1, k - 1] + Fs_M[i - 1, j - 1, k
           - 1])/2 else 0;
        Ft_y[i - 1, j - 1, k - 1] = if convective then (dy[j -
          1]*Ft_M[i - 1, j - 1, k - 1] + dy[j]*Ft_M[i - 1, j,
          k - 1])/(dy[j - 1] + dy[j]) else 0;
        Fb_y[i - 1, j - 1, k - 1] = if convective then (dy[j -
          1]*Fb_M[i - 1, j - 1, k - 1] + dy[j]*Fb_M[i - 1, j,
          k - 1])/(dy[j - 1] + dy[j]) else 0;

        // Peclet number
        Pe_y[i - 1, j - 1, k - 1] = Fe_y[i - 1, j - 1,
          k - 1]/De_y[i - 1, j - 1, k - 1];
        Pw_y[i - 1, j - 1, k - 1] = Fw_y[i - 1, j - 1,
          k - 1]/Dw_y[i - 1, j - 1, k - 1];
        Pn_y[i - 1, j - 1, k - 1] = Fn_y[i - 1, j - 1,
          k - 1]/Dn_y[i - 1, j - 1, k - 1];
        Ps_y[i - 1, j - 1, k - 1] = Fs_y[i - 1, j - 1,
          k - 1]/Ds_y[i - 1, j - 1, k - 1];
        Pt_y[i - 1, j - 1, k - 1] = Ft_y[i - 1, j - 1,
          k - 1]/Dt_y[i - 1, j - 1, k - 1];
        Pb_y[i - 1, j - 1, k - 1] = Fb_y[i - 1, j - 1,
          k - 1]/Db_y[i - 1, j - 1, k - 1];

        // terms of the momentum equation
        aE_y[i - 1, j - 1, k - 1] = De_y[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pe_y[i - 1, j - 1, k - 1], intScheme) + max(-
          Fe_y[i - 1, j - 1, k - 1], 0);
        aW_y[i - 1, j - 1, k - 1] = Dw_y[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pw_y[i - 1, j - 1, k - 1], intScheme) + max(
          Fw_y[i - 1, j - 1, k - 1], 0);
        aN_y[i - 1, j - 1, k - 1] = Dn_y[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pn_y[i - 1, j - 1, k - 1], intScheme) + max(-
          Fn_y[i - 1, j - 1, k - 1], 0);
        aS_y[i - 1, j - 1, k - 1] = Ds_y[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Ps_y[i - 1, j - 1, k - 1], intScheme) + max(
          Fs_y[i - 1, j - 1, k - 1], 0);
        aT_y[i - 1, j - 1, k - 1] = Dt_y[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pt_y[i - 1, j - 1, k - 1], intScheme) + max(-
          Ft_y[i - 1, j - 1, k - 1], 0);
        aB_y[i - 1, j - 1, k - 1] = Db_y[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pb_y[i - 1, j - 1, k - 1], intScheme) + max(
          Fb_y[i - 1, j - 1, k - 1], 0);

        // Here the terms
        // + Fe_y[i,J,k] - Fw_y[i,J,k] + Fn_y[i,J,k] - Fs_y[i,J,k] + Ft_y[i,J,k] - Fb_y[i,J,k]
        // to be added are disregarded since the mass variation is not taken into account
        aP_y[i - 1, j - 1, k - 1] = aE_y[i - 1, j - 1, k - 1]
           + aW_y[i - 1, j - 1, k - 1] + aN_y[i - 1, j - 1,
          k - 1] + aS_y[i - 1, j - 1, k -
          1] + aB_y[i - 1, j - 1, k - 1] + aT_y[i - 1, j - 1,
          k - 1];

        // y-momentum equation
        (dy[j - 1]*0.5 + dy[j]*0.5)*dx[i - 1]*dz[k - 1]*
          rho_o*der(Vy[i, j, k]) + aP_y[i - 1, j - 1,
          k - 1]*Vy[i, j, k] = aE_y[i - 1,
          j - 1, k - 1]*Vy[i + 1, j, k] +
          aW_y[i - 1, j - 1, k - 1]*Vy[i - 1, j, k]
           + aN_y[i - 1, j - 1, k - 1]*Vy[i, j + 1, k]
           + aS_y[i - 1, j - 1, k - 1]*Vy[i, j - 1, k]
           + aT_y[i - 1, j - 1, k - 1]*Vy[i, j, k
           + 1] + aB_y[i - 1, j - 1, k - 1]*Vy[i, j,
          k - 1] + (P[i - 1, j - 1, k - 1]
           - P[i - 1, j, k - 1])*dx[i - 1]*dz[k
           - 1];

      end for;
    end for;
  end for;
  // END Y-MOMENTUM PRESERVATION
  ////////////////////////////////////////

  ////////////////////////////////////////
  // Z-MOMENTUM PRESERVATION
  for i in 2:I+1 loop
    for j in 2:J+1 loop
      for k in 2:K loop

        // diffusive fluxes
        De_z[i - 1, j - 1, k - 1] = dy[j - 1]*(dz[k
           - 1]*0.5 + dz[k]*0.5)*(0.5*MUz[i, j, k]
           + 0.5*MUz[i + 1, j, k])/dxz[i];
        Dw_z[i - 1, j - 1, k - 1] = dy[j - 1]*(dz[k
           - 1]*0.5 + dz[k]*0.5)*(0.5*MUz[i, j, k]
           + 0.5*MUz[i - 1, j, k])/dxz[i - 1];
        Dn_z[i - 1, j - 1, k - 1] = dx[i - 1]*(dz[k
           - 1]*0.5 + dz[k]*0.5)*(0.5*MUz[i, j, k]
           + 0.5*MUz[i, j + 1, k])/dyz[j];
        Ds_z[i - 1, j - 1, k - 1] = dx[i - 1]*(dz[k
           - 1]*0.5 + dz[k]*0.5)*(0.5*MUz[i, j, k]
           + 0.5*MUz[i, j - 1, k])/dyz[j - 1];
        Dt_z[i - 1, j - 1, k - 1] = dy[j - 1]*dx[i - 1]*(0.5*
          MUz[i, j, k] + 0.5*MUz[i, j, k +
          1])/dz[k];
        Db_z[i - 1, j - 1, k - 1] = dy[j - 1]*dx[i - 1]*(0.5*
          MUz[i, j, k] + 0.5*MUz[i, j, k -
          1])/dz[k - 1];

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fe_z[i - 1, j - 1, k - 1] = if convective then (dz[k
           - 1]*Fe_M[i - 1, j - 1, k - 1] + dz[k]
          *Fe_M[i - 1, j - 1, k])/(dz[k -
          1] + dz[k]) else 0;
        Fw_z[i - 1, j - 1, k - 1] = if convective then (dz[
          k - 1]*Fw_M[i - 1, j - 1, k - 1]
           + dz[k]*Fw_M[i - 1, j - 1, k])/
          (dz[k - 1] + dz[k]) else 0;
        Fn_z[i - 1, j - 1, k - 1] = if convective then (dz[
          k - 1]*Fe_M[i - 1, j - 1, k - 1]
           + dz[k]*Fe_M[i - 1, j - 1, k])/
          (dz[k - 1] + dz[k]) else 0;
        Fs_z[i - 1, j - 1, k - 1] = if convective then (dz[
          k - 1]*Fw_M[i - 1, j - 1, k - 1]
           + dz[k]*Fw_M[i - 1, j - 1, k])/
          (dz[k - 1] + dz[k]) else 0;
        Ft_z[i - 1, j - 1, k - 1] = if convective then (Fb_M[
          i - 1, j - 1, k] + Ft_M[i - 1, j - 1, k])
          /2 else 0;
        Fb_z[i - 1, j - 1, k - 1] = if convective then (Fb_M[
          i - 1, j - 1, k - 1] + Ft_M[i - 1, j - 1, k
           - 1])/2 else 0;

        // Peclet number
        Pe_z[i - 1, j - 1, k - 1] = Fe_z[i - 1, j - 1,
          k - 1]/De_z[i - 1, j - 1, k - 1];
        Pw_z[i - 1, j - 1, k - 1] = Fw_z[i - 1, j - 1,
          k - 1]/Dw_z[i - 1, j - 1, k - 1];
        Pn_z[i - 1, j - 1, k - 1] = Fn_z[i - 1, j - 1,
          k - 1]/Dn_z[i - 1, j - 1, k - 1];
        Ps_z[i - 1, j - 1, k - 1] = Fs_z[i - 1, j - 1,
          k - 1]/Ds_z[i - 1, j - 1, k - 1];
        Pt_z[i - 1, j - 1, k - 1] = Ft_z[i - 1, j - 1,
          k - 1]/Dt_z[i - 1, j - 1, k - 1];
        Pb_z[i - 1, j - 1, k - 1] = Fb_z[i - 1, j - 1,
          k - 1]/Db_z[i - 1, j - 1, k - 1];

        // terms of the momentum equation
        aE_z[i - 1, j - 1, k - 1] = De_z[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pe_z[i - 1, j - 1, k - 1], intScheme) + max(-
          Fe_z[i - 1, j - 1, k - 1], 0);
        aW_z[i - 1, j - 1, k - 1] = Dw_z[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pw_z[i - 1, j - 1, k - 1], intScheme) + max(
          Fw_z[i - 1, j - 1, k - 1], 0);
        aN_z[i - 1, j - 1, k - 1] = Dn_z[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pn_z[i - 1, j - 1, k - 1], intScheme) + max(-
          Fn_z[i - 1, j - 1, k - 1], 0);
        aS_z[i - 1, j - 1, k - 1] = Ds_z[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Ps_z[i - 1, j - 1, k - 1], intScheme) + max(
          Fs_z[i - 1, j - 1, k - 1], 0);
        aT_z[i - 1, j - 1, k - 1] = Dt_z[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pt_z[i - 1, j - 1, k - 1], intScheme) + max(-
          Ft_z[i - 1, j - 1, k - 1], 0);
        aB_z[i - 1, j - 1, k - 1] = Db_z[i - 1, j - 1, k - 1]*
          CFDModelica.Functions.A(Pb_z[i - 1, j - 1, k - 1], intScheme) + max(
          Fb_z[i - 1, j - 1, k - 1], 0);

        // Here the terms
        // Fe_z[i,J,k] - Fw_z[i,J,k] + Fn_z[i,J,k] - Fs_z[i,J,k] + Ft_z[i,J,k] - Fb_z[i,J,k]
        // to be added are disregarded since the mass variation is not taken into account
        aP_z[i - 1, j - 1, k - 1] = aE_z[i - 1, j - 1, k - 1]
           + aW_z[i - 1, j - 1, k - 1] + aN_z[i - 1, j - 1,
          k - 1] + aS_z[i - 1, j - 1, k -
          1] + aB_z[i - 1, j - 1, k - 1] + aT_z[i - 1, j - 1,
          k - 1];

        // z-momentum equation
        (dz[k - 1]*0.5 + dz[k]*0.5)*dy[j -
          1]*dx[i - 1]*rho_o*der(Vz[i, j, k]) + aP_z[i - 1,
          j - 1, k - 1]*Vz[i, j, k] =
          aE_z[i - 1, j - 1, k - 1]*Vz[i + 1, j, k]
           + aW_z[i - 1, j - 1, k - 1]*Vz[i - 1, j, k]
           + aN_z[i - 1, j - 1, k - 1]*Vz[i, j + 1, k]
           + aS_z[i - 1, j - 1, k - 1]*Vz[i, j - 1, k]
           + aT_z[i - 1, j - 1, k - 1]*Vz[i, j, k
           + 1] + aB_z[i - 1, j - 1, k - 1]*Vz[i, j,
          k - 1] + (P[i - 1, j - 1, k - 1]
           - P[i - 1, j - 1, k])*dy[j - 1]*dx[i - 1] - (if
          gravity then (rho[i, j, k] + rho[i, j, k
           + 1])/2*(dz[k - 1]*0.5 + dz[k]*
          0.5)*dy[j - 1]*dx[i - 1]*g else 0);

      end for;
    end for;
  end for;
  // END Z-MOMENTUM PRESERVATION
  ////////////////////////////////////////
end Room3D_CFDwall_momentum;
