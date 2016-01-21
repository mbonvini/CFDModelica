within CFDModelica.Flow3D;
partial model Room3D_CFDwall_momentum_boundaries
  extends Room3D_CFDwall_base;
equation
  ////////////////////////////////////////
  // X-MOMENTUM PRESERVATION
  // LEFT-WALL
  for j in 2:J+1 loop
    for k in 2:K + 1 loop

      if CFDModelica.Functions.isIn(
          j - 1,
          k - 1,
          jk_left,
          all_left) then
        // diffusive fluxes
        De_bl_x[j - 1, k - 1] = dy[j - 1]*dz[k
           - 1]*(0.5*MUx[1, j, k] + 0.5*MUx[2, j, k])
          /(dx[1]/2);

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fe_bl_x[j - 1, k - 1] = if convective then (Fe_M[1, j - 1,
          k - 1] + Fw_M[1, j - 1, k - 1])/
          2 else 0;

        // Peclet number
        Pe_bl_x[j - 1, k - 1] = Fe_bl_x[j - 1, k
           - 1]/De_bl_x[j - 1, k - 1];

        // terms of the momentum equation
        aE_bl_x[j - 1, k - 1] = De_bl_x[j - 1, k - 1]*CFDModelica.Functions.A(
          Pe_bl_x[j - 1, k - 1], intScheme) + max(-Fe_bl_x[j - 1, k - 1], 0);

        // x-momentum equation
        /*
        (dx[1]*0.5)*dy[j-1]*dz[k-1]*rho_o*der(Vx[1,j,k]) +
        aE_bl_x[j-1,k-1]*Vx[1,j,k]
        =
        aE_bl_x[j-1,k-1]*Vx[2,j,k] +
        (Pleft[j-1,k-1] - P[1,j-1,k-1])*dy[j-1]*dz[k-1];
        */
        Pleft[j - 1, k - 1] = P[1, j - 1, k
           - 1];
      else
        De_bl_x[j - 1, k - 1] = 0;
        Fe_bl_x[j - 1, k - 1] = 0;
        Pe_bl_x[j - 1, k - 1] = 0;
        aE_bl_x[j - 1, k - 1] = 0;
        Pleft[j - 1, k - 1] = 0;
      end if;

    end for;
  end for;
  // RIGHT-WALL
  for j in 2:J+1 loop
    for k in 2:K + 1 loop

      if CFDModelica.Functions.isIn(
          j - 1,
          k - 1,
          jk_right,
          all_right) then
      // diffusive fluxes
        Dw_br_x[j - 1, k - 1] = dy[j - 1]*dz[k
           - 1]*(0.5*MUx[I, j, k] + 0.5*MUx[I + 1, j,
          k])/dx[I]/2;

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fw_br_x[j - 1, k - 1] = if convective then (Fw_M[I, j - 1,
          k - 1] + Fe_M[I, j - 1, k - 1])/
          2 else 0;

        // Peclet number
        Pw_br_x[j - 1, k - 1] = Fw_br_x[j - 1, k
           - 1]/Dw_br_x[j - 1, k - 1];

        // terms of the momentum equation
        aW_br_x[j - 1, k - 1] = Dw_br_x[j - 1, k - 1]*CFDModelica.Functions.A(
          Pw_br_x[j - 1, k - 1], intScheme) + max(Fw_br_x[j - 1, k - 1], 0);

        // x-momentum equation
        /*
        (dx[I]*0.5)*dy[j-1]*dz[k-1]*rho_o*der(Vx[I+1,j,k]) +
        aW_br_x[j-1,k-1]*Vx[I+1,j,k]
        =
        aW_br_x[j-1,k-1]*Vx[I,j,k] +
        (P[I,j-1,k-1] - Pright[j-1,k-1])*dy[j-1]*dz[k-1];
        */
        P[I, j - 1, k - 1] = Pright[j - 1, k
           - 1];
      else
        Dw_br_x[j - 1, k - 1] = 0;
        Fw_br_x[j - 1, k - 1] = 0;
        Pw_br_x[j - 1, k - 1] = 0;
        aW_br_x[j - 1, k - 1] = 0;
        Pright[j - 1, k - 1] = 0;
      end if;

    end for;
  end for;
  // END X-MOMENTUM PRESERVATION
  ////////////////////////////////////////

  ////////////////////////////////////////
  // Y-MOMENTUM PRESERVATION
  // FRONT
  for i in 2:I+1 loop
    for k in 2:K + 1 loop

      if CFDModelica.Functions.isIn(
          i - 1,
          k - 1,
          ik_front,
          all_front) then
      // diffusive fluxes
        Dn_bf_y[i - 1, k - 1] = dx[i - 1]*dz[k
           - 1]*(0.5*MUy[i, 1, k] + 0.5*MUy[i, 2, k])
          /dy[1]/2;

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fn_bf_y[i - 1, k - 1] = if convective then (Fn_M[i - 1, 1,
          k - 1] + Fs_M[i - 1, 1, k - 1])/
          2 else 0;

        // Peclet number
        Pn_bf_y[i - 1, k - 1] = Fn_bf_y[i - 1, k
           - 1]/Dn_bf_y[i - 1, k - 1];

        // terms of the momentum equation
        aN_bf_y[i - 1, k - 1] = Dn_bf_y[i - 1, k - 1]*CFDModelica.Functions.A(
          Pn_bf_y[i - 1, k - 1], intScheme) + max(-Fn_bf_y[i - 1, k - 1], 0);

        /*
        // y-momentum equation
        (dy[1]*0.5)*dx[i-1]*dz[k-1]*rho_o*der(Vy[i,1,k]) +
        aN_bf_y[i-1,k-1]*Vy[i,1,k]
        =
        aN_bf_y[i-1,k-1]*Vy[i,2,k] +
        (Pfront[i-1,k-1] - P[i-1,1,k-1])*dx[i-1]*dz[k-1];
        */
        P[i - 1, 1, k - 1] = Pfront[i - 1, k
           - 1];
      else
        Dn_bf_y[i - 1, k - 1] = 0;
        Fn_bf_y[i - 1, k - 1] = 0;
        Pn_bf_y[i - 1, k - 1] = 0;
        aN_bf_y[i - 1, k - 1] = 0;
        Pfront[i - 1, k - 1] = 0;
      end if;

    end for;
  end for;
  //REAR
  for i in 2:I+1 loop
    for k in 2:K + 1 loop

      if CFDModelica.Functions.isIn(
          i - 1,
          k - 1,
          ik_rear,
          all_rear) then
      // diffusive fluxes
        Ds_bre_y[i - 1, k - 1] = dx[i - 1]*dz[k
           - 1]*(0.5*MUy[i, J, k] + 0.5*MUy[i, J + 1,
          k])/dy[J]/2;

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fs_bre_y[i - 1, k - 1] = if convective then (Fs_M[i - 1, J,
          k - 1] + Fn_M[i - 1, J, k - 1])/
          2 else 0;

        // Peclet number
        Ps_bre_y[i - 1, k - 1] = Fs_bre_y[i - 1, k
           - 1]/Ds_bre_y[i - 1, k - 1];

        // terms of the momentum equation
        aS_bre_y[i - 1, k - 1] = Ds_bre_y[i - 1, k - 1]*CFDModelica.Functions.A(
          Ps_bre_y[i - 1, k - 1], intScheme) + max(Fs_bre_y[i - 1, k - 1], 0);

        /*
        // y-momentum equation
        (dy[J]*0.5)*dx[i-1]*dz[k-1]*rho_o*der(Vy[i,J+1,k]) +
        aS_bre_y[i-1,k-1]*Vy[i,J+1,k]
        =
        aS_bre_y[i-1,k-1]*Vy[i,J,k] +
        (P[i-1,J,k-1] - Prear[i-1,k-1])*dx[i-1]*dz[k-1];
        */
        P[i - 1, J, k - 1] = Prear[i - 1, k
           - 1];
      else
        Ds_bre_y[i - 1, k - 1] = 0;
        Fs_bre_y[i - 1, k - 1] = 0;
        Ps_bre_y[i - 1, k - 1] = 0;
        aS_bre_y[i - 1, k - 1] = 0;
        Prear[i - 1, k - 1] = 0;
      end if;

    end for;
  end for;
  // END Y-MOMENTUM PRESERVATION
  ////////////////////////////////////////

  ////////////////////////////////////////
  // Z-MOMENTUM PRESERVATION
  // BOTTOM
  for i in 2:I+1 loop
    for j in 2:J+1 loop

      if CFDModelica.Functions.isIn(
          i - 1,
          j - 1,
          ij_floor,
          all_floor) then
      // diffusive fluxes
        Dt_bb_z[i-1,j-1] = dy[j-1]*dx[i-1]*(0.5*MUz[i,j,1]+0.5*MUz[i,j,2])/dz[1]/2;

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Ft_bb_z[i-1,j-1] = if convective then (Ft_M[i-1,j-1,1]+Fb_M[i-1,j-1,1])/2 else 0;

        // Peclet number
        Pt_bb_z[i-1,j-1] = Ft_bb_z[i-1,j-1]/Dt_bb_z[i-1,j-1];

        // terms of the momentum equation
        aT_bb_z[i - 1, j - 1] = Dt_bb_z[i - 1, j - 1]*CFDModelica.Functions.A(
          Pt_bb_z[i - 1, j - 1], intScheme) + max(-Ft_bb_z[i - 1, j - 1], 0);

        /*
        // z-momentum equation
        (dz[1]*0.5)*dy[j-1]*dx[i-1]*rho_o*der(Vz[i,j,1]) +
        aT_bb_z[i-1,j-1]*Vz[i,j,1]
        =
        aT_bb_z[i-1,j-1]*Vz[i,j,2] +
        (Pfloor[i-1,j-1] - P[i-1,j-1,1])*dy[j-1]*dx[i-1] -
        (if gravity then (rho[i,j,1]+rho[i,j,2])/2*(dz[1]*0.5)*dy[j-1]*dx[i-1]*g else 0);
        */
        Pfloor[i-1,j-1] = P[i-1,j-1,1];
      else
        Dt_bb_z[i-1,j-1] = 0;
        Ft_bb_z[i-1,j-1] = 0;
        Pt_bb_z[i-1,j-1] = 0;
        aT_bb_z[i-1,j-1] = 0;
        Pfloor[i-1,j-1]  = 0;
      end if;

    end for;
  end for;
  // TOP
  for i in 2:I+1 loop
    for j in 2:J+1 loop

      if CFDModelica.Functions.isIn(
          i - 1,
          j - 1,
          ij_cei,
          all_cei) then
      // diffusive fluxes
        Db_bt_z[i-1,j-1] = dy[j-1]*dx[i-1]*(0.5*MUz[i,j,K]+0.5*MUz[i,j,K+1])/dz[K]/2;

        // mass fluxes
        // is exactly the algebraic mean, since the faces are located in the middle
        Fb_bt_z[i-1,j-1] = if convective then (Fb_M[i-1,j-1,K]+Ft_M[i-1,j-1,K])/2 else 0;

        // Peclet number
        Pb_bt_z[i-1,j-1] = Fb_bt_z[i-1,j-1]/Db_bt_z[i-1,j-1];

        // terms of the momentum equation
        aB_bt_z[i - 1, j - 1] = Db_bt_z[i - 1, j - 1]*CFDModelica.Functions.A(
          Pb_bt_z[i - 1, j - 1], intScheme) + max(Fb_bt_z[i - 1, j - 1], 0);

        // z-momentum equation
        /*
        (dz[K]*0.5)*dy[j-1]*dx[i-1]*rho_o*der(Vz[i,j,K+1]) +
        aB_bt_z[i-1,j-1]*Vz[i,j,K+1]
        =
        aB_bt_z[i-1,j-1]*Vz[i,j,K] +
        (P[i-1,j-1,K] - Pcei[i-1,j-1])*dy[j-1]*dx[i-1] -
        (if gravity then (rho[i,j,K+1]+rho[i,j,K+2])/2*(dz[K]*0.5)*dy[j-1]*dx[i-1]*g else 0);
        */
        P[i-1,j-1,K] = Pcei[i-1,j-1];
      else
        Db_bt_z[i-1,j-1] = 0;
        Fb_bt_z[i-1,j-1] = 0;
        Pb_bt_z[i-1,j-1] = 0;
        aB_bt_z[i-1,j-1] = 0;
        Pcei[i-1,j-1]    = 0;
      end if;

    end for;
  end for;
  // END Z-MOMENTUM PRESERVATION
  ////////////////////////////////////////
end Room3D_CFDwall_momentum_boundaries;
