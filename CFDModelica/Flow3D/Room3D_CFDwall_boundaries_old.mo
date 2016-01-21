within CFDModelica.Flow3D;
partial model Room3D_CFDwall_boundaries_old
  extends Room3D_CFDwall_base;
equation
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // BOUNDARY VALUES FOR:
  // X-VELOCITIES
  for k in 1:K loop
    for j in 1:J loop

      // left boundaries (vx)
      Vx[1, j + 1, k + 1] = +fluidPort_left[j, k].m_flow/(dy[
        j]*dz[k])/rho_o;                            //rho[1,j+1,k+1];
      // right boundaries (vx)
      Vx[I + 1, j + 1, k + 1] = -fluidPort_right[j, k].m_flow
        /(dy[j]*dz[k])/rho_o;                           //rho[I+2,j+1,k+1];

    end for;
  end for;

  for i in 1:I+1 loop
    for j in 1:J+2 loop
      // top and bottom boundaries (vx)
      Vx[i,j,1]    = 0;
      Vx[i,j,K+2]  = 0;
    end for;
  end for;

  for i in 1:I+1 loop
    for k in 2:K + 1 loop
      // front and rear boundaries (vx)
      Vx[i, 1, k] = 0;
      Vx[i, J + 2, k] = 0;
    end for;
  end for;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // BOUNDARY VALUES FOR:
  // Y-VELOCITIES
  for i in 1:I loop
    for k in 1:K loop

      // front boundaries (vy)
      Vy[i + 1, 1, k + 1] = +fluidPort_front[i, k].m_flow/(
        dx[i]*dz[k])/rho_o;                          //rho[i+1,1,k+1];
      // rear boundaries (vy)
      Vy[i + 1, J + 1, k + 1] = -fluidPort_rear[i, k].m_flow/
        (dx[i]*dz[k])/rho_o;                          //rho[i+1,J+2,k+1];

    end for;
  end for;

  for i in 1:I+2 loop
    for j in 1:J+1 loop
      // bottom and top (vy)
      Vy[i,j,1]   = 0;
      Vy[i,j,K+2] = 0;
    end for;
  end for;

  for j in 1:J+1 loop
    for k in 2:K + 1 loop
      // left and right (vy)
      Vy[1, j, k] = 0;
      Vy[I + 2, j, k] = 0;
    end for;
  end for;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // BOUNDARY VALUES FOR:
  // Z-VELOCITIES
  for i in 1:I loop
    for j in 1:J loop

      // bottom boundaries (vz)
      Vz[i+1,j+1,1] = + fluidPort_floor[i,j].m_flow/(dx[i]*dy[j])/rho_o;//rho[i+1,j+1,1];
      // top boundaries (vz)
      Vz[i+1,j+1,K+1] = - fluidPort_cei[i,j].m_flow/(dx[i]*dy[j])/rho_o;//rho[i+1,j+1,K+2];

    end for;
  end for;

  for k in 1:K + 1 loop
    for j in 1:J+2 loop
      // left and right (vz)
      Vz[1, j, k] = 0;
      Vz[I + 2, j, k] = 0;
    end for;
  end for;

  for i in 2:I+1 loop
    for k in 1:K + 1 loop
      // front and rear (vz)
      Vz[i, 1, k] = 0;
      Vz[i, J + 2, k] = 0;
    end for;
  end for;

  ///////////////////////////////////////////////////////
  // BOUNDARY CONDITIONS FOR TEMPERATURES AND DENSITIES
  // the eight corners are not relevant

  // (1,1,1) corner
  T[1,1,1]         = 0;
  rho[1,1,1]       = 0;

  // (I+2,1,1) corner
  T[I+2,1,1]       = 0;
  rho[I+2,1,1]     = 0;

  // (1,1,K+2) corner
  T[1,1,K+2]       = 0;
  rho[1,1,K+2]     = 0;

  // (I+2,1,K+2) corner
  T[I+2,1,K+2]     = 0;
  rho[I+2,1,K+2]   = 0;

  // (1,J+2,1) corner
  T[1,J+2,1]       = 0;
  rho[1,J+2,1]     = 0;

  // (I+2,J+2,1) corner
  T[I+2,J+2,1]     = 0;
  rho[I+2,J+2,1]   = 0;

  // (1,J+2,K+2) corner
  T[1,J+2,K+2]     = 0;
  rho[1,J+2,K+2]   = 0;

  // (I+2,J+2,K+2) corner
  T[I+2,J+2,K+2]   = 0;
  rho[I+2,J+2,K+2] = 0;

  // the twelve edges are not relevant
  for i in 1:I loop
    T[i+1,1,1]       = 0;
    rho[i+1,1,1]     = 0;
    T[i+1,J+2,1]     = 0;
    rho[i+1,J+2,1]   = 0;
    T[i+1,1,K+2]     = 0;
    rho[i+1,1,K+2]   = 0;
    T[i+1,J+2,K+2]   = 0;
    rho[i+1,J+2,K+2] = 0;
  end for;
  for j in 1:J loop
    T[1,j+1,1]       = 0;
    rho[1,j+1,1]     = 0;
    T[I+2,j+1,1]     = 0;
    rho[I+2,j+1,1]   = 0;
    T[1,j+1,K+2]     = 0;
    rho[1,j+1,K+2]   = 0;
    T[I+2,j+1,K+2]   = 0;
    rho[I+2,j+1,K+2] = 0;
  end for;
  for k in 1:K loop
    T[1, 1, k + 1] = 0;
    rho[1, 1, k + 1] = 0;
    T[I + 2, 1, k + 1] = 0;
    rho[I + 2, 1, k + 1] = 0;
    T[1, J + 2, k + 1] = 0;
    rho[1, J + 2, k + 1] = 0;
    T[I + 2, J + 2, k + 1] = 0;
    rho[I + 2, J + 2, k + 1] = 0;
  end for;

  /******************************************************************/
  // the wall surfaces
  /******************************************************************/
  // left and right
  for j in 2:J+1 loop
    for k in 2:K + 1 loop

      // densities, enthalpies and boundary pressures
      // left
      if CFDModelica.Functions.isIn(
          j - 1,
          k - 1,
          jk_left,
          all_left) then
        rho[1, j, k] = Medium.density_ph(Pleft[j - 1,
          k - 1], actualStream(fluidPort_left[j - 1, k - 1].h_outflow));
        fluidPort_left[j - 1, k - 1].h_outflow = Medium.specificEnthalpy_pT(P[1,
          j - 1, k - 1], T[2, j, k]);
        Pleft[j - 1, k - 1] = fluidPort_left[j - 1, k - 1].p -
          Medium.p_ref;
      else
        rho[1, j, k] = rho_o;
        fluidPort_left[j-1,k-1].h_outflow  = 0;
        fluidPort_left[j-1,k-1].p = Medium.p_ref;
      end if;
      // densities, enthalpies and boundary pressures
      // right
      if CFDModelica.Functions.isIn(
          j - 1,
          k - 1,
          jk_right,
          all_right) then
        rho[I + 2, j, k] = Medium.density_ph(Pright[j - 1,
          k - 1], actualStream(fluidPort_right[j - 1, k - 1].h_outflow));
        fluidPort_right[j - 1, k - 1].h_outflow = Medium.specificEnthalpy_pT(P[
          I, j - 1, k - 1], T[I + 1, j, k]);
        Pright[j - 1, k - 1] = fluidPort_right[j - 1, k - 1].p
           - Medium.p_ref;
      else
        rho[I + 2, j, k] = rho_o;
        fluidPort_right[j-1,k-1].h_outflow = 0;
        fluidPort_right[j-1,k-1].p = Medium.p_ref;
      end if;

      // Mixing temperatures
      T[1, j, k] = (alpha*abs(heatPort_left[j - 1, k - 1].Q_flow)
        *heatPort_left[j - 1, k - 1].T + beta*max(0, fluidPort_left[j - 1, k - 1].m_flow)
        *Medium.temperature_ph(Pleft[j - 1, k - 1], inStream(
         fluidPort_left[j - 1, k - 1].h_outflow)) + Modelica.Constants.eps*T[2,
        j, k])/(alpha*abs(heatPort_left[j - 1, k - 1].Q_flow)
         + beta*max(0, fluidPort_left[j - 1, k - 1].m_flow) + Modelica.Constants.eps);                                                                                                    //Tmix_left[j-1,k-1];
      T[I + 2, j, k] = (alpha*abs(heatPort_right[j - 1, k - 1].Q_flow)*heatPort_right[j - 1, k - 1].T + beta*max(0, fluidPort_right[j - 1, k - 1].m_flow)*Medium.temperature_ph(Pright[j - 1, k - 1], inStream(fluidPort_right[j - 1, k - 1].h_outflow)) + Modelica.Constants.eps*T[I + 1, j, k])/(alpha*abs(
        heatPort_right[j - 1, k - 1].Q_flow) + beta*max(0, fluidPort_right[j - 1,
        k - 1].m_flow) + Modelica.Constants.eps);                                                                                                    //Tmix_right[j-1,k-1];

      // heat flows
      heatPort_left[j - 1, k - 1].Q_flow = Dw_M[1, j - 1, k -
        1]*(heatPort_left[j - 1, k - 1].T - T[2, j, k]);
      heatPort_right[j - 1, k - 1].Q_flow = De_M[I, j - 1, k -
        1]*(heatPort_right[j - 1, k - 1].T - T[I + 1, j, k]);

      end for;
  end for;

  // front and rear
  for i in 2:I+1 loop
    for k in 2:K + 1 loop

      // densities, enthalpies and boundary pressures
      // front
      if CFDModelica.Functions.isIn(
          i - 1,
          k - 1,
          ik_front,
          all_front) then
        rho[i, 1, k] = Medium.density_ph(Pfront[i - 1,
          k - 1], actualStream(fluidPort_front[i - 1, k - 1].h_outflow));
        fluidPort_front[i - 1, k - 1].h_outflow = Medium.specificEnthalpy_pT(P[
          i - 1, 1, k - 1], T[i, 2, k]);
        Pfront[i - 1, k - 1] = fluidPort_front[i - 1, k - 1].p
           - Medium.p_ref;
      else
        rho[i, 1, k] = rho_o;
        fluidPort_front[i-1,k-1].h_outflow = 0;
        fluidPort_front[i-1,k-1].p = Medium.p_ref;
      end if;
      // rear
      if CFDModelica.Functions.isIn(
          i - 1,
          k - 1,
          ik_rear,
          all_rear) then
        rho[i, J + 2, k] = Medium.density_ph(Prear[i - 1,
          k - 1], actualStream(fluidPort_rear[i - 1, k - 1].h_outflow));
        fluidPort_rear[i - 1, k - 1].h_outflow = Medium.specificEnthalpy_pT(P[i -
          1, J, k - 1], T[i, J + 1, k]);
        Prear[i - 1, k - 1] = fluidPort_rear[i - 1, k - 1].p -
          Medium.p_ref;
      else
        rho[i, J + 2, k] = rho_o;
        fluidPort_rear[i-1,k-1].h_outflow  = 0;
        fluidPort_rear[i-1,k-1].p = Medium.p_ref;
      end if;

      // Mixing temperatures
      T[i, 1, k] = (alpha*abs(heatPort_front[i - 1, k - 1].Q_flow)
        *heatPort_front[i - 1, k - 1].T + beta*max(0, fluidPort_front[i - 1, k -
        1].m_flow)*Medium.temperature_ph(Pfront[i - 1, k - 1],
        inStream(fluidPort_front[i - 1, k - 1].h_outflow)) + Modelica.Constants.eps
        *T[i, 2, k])/(alpha*abs(heatPort_front[i - 1, k - 1].Q_flow)
         + beta*max(0, fluidPort_front[i - 1, k - 1].m_flow) + Modelica.Constants.eps);                                                                                                    //Tmix_front[i-1,k-1];
      T[i, J + 2, k] = (alpha*abs(heatPort_rear[i - 1, k - 1].Q_flow)*heatPort_rear[i - 1, k - 1].T + beta*max(0, fluidPort_rear[i - 1, k - 1].m_flow)*Medium.temperature_ph(Prear[i - 1, k - 1], inStream(fluidPort_rear[i - 1, k - 1].h_outflow)) + Modelica.Constants.eps*T[i, J + 1, k])/(alpha*abs(
        heatPort_rear[i - 1, k - 1].Q_flow) + beta*max(0, fluidPort_rear[i - 1,
        k - 1].m_flow) + Modelica.Constants.eps);                                                                                                    //Tmix_rear[i-1,k-1];

      // heat flows
      heatPort_front[i - 1, k - 1].Q_flow = Ds_M[i - 1, 1, k -
        1]*(heatPort_front[i - 1, k - 1].T - T[i, 2, k]);
      heatPort_rear[i - 1, k - 1].Q_flow = Dn_M[i - 1, J, k -
        1]*(heatPort_rear[i - 1, k - 1].T - T[i, J + 1, k]);

      end for;
  end for;

  // floor and ceiling
  for i in 2:I+1 loop
    for j in 2:J+1 loop

      // densities, enthalpies and boundary pressures
      // bottom
      if CFDModelica.Functions.isIn(
          i - 1,
          j - 1,
          ij_floor,
          all_floor) then
        rho[i,j,1] = Medium.density_ph(Pfloor[i-1,j-1],actualStream(fluidPort_floor[i-1,j-1].h_outflow));
        fluidPort_floor[i-1,j-1].h_outflow = Medium.specificEnthalpy_pT(P[i-1,j-1,1],T[i,j,2]);
        Pfloor[i-1,j-1] = fluidPort_floor[i-1,j-1].p - Medium.p_ref;
      else
        rho[i,j,1] = rho_o;
        fluidPort_floor[i-1,j-1].h_outflow = 0;
        fluidPort_floor[i-1,j-1].p = Medium.p_ref;
      end if;
      // top
      if CFDModelica.Functions.isIn(
          i - 1,
          j - 1,
          ij_cei,
          all_cei) then
        rho[i,j,K+2] = Medium.density_ph(Pcei[i-1,j-1],actualStream(fluidPort_cei[i-1,j-1].h_outflow));
        fluidPort_cei[i-1,j-1].h_outflow = Medium.specificEnthalpy_pT(P[i-1,j-1,K],T[i,j,K+1]);
        Pcei[i-1,j-1]   = fluidPort_cei[i-1,j-1].p - Medium.p_ref;
      else
        rho[i,j,K+2] = rho_o;
        fluidPort_cei[i-1,j-1].h_outflow = 0;
        fluidPort_cei[i-1,j-1].p = Medium.p_ref;
      end if;

      // Mixing temperatures
      T[i, j, 1] = (alpha*abs(heatPort_floor[i - 1, j - 1].Q_flow)*
        heatPort_floor[i - 1, j - 1].T + beta*max(0, fluidPort_floor[i - 1, j -
        1].m_flow)*Medium.temperature_ph(Pfloor[i - 1, j - 1], inStream(
        fluidPort_floor[i - 1, j - 1].h_outflow)) + Modelica.Constants.eps*T[i,
        j, 2])/(alpha*abs(heatPort_floor[i - 1, j - 1].Q_flow) + beta*max(0,
        fluidPort_floor[i - 1, j - 1].m_flow) + Modelica.Constants.eps);                                                                                                    //Tmix_cei[i-1,j-1];
      T[i, j, K + 2] = (alpha*abs(heatPort_cei[i - 1, j - 1].Q_flow)*heatPort_cei[i - 1, j - 1].T + beta*max(0, fluidPort_cei[i - 1, j - 1].m_flow)*Medium.temperature_ph(Pcei[i - 1, j - 1], inStream(fluidPort_cei[i - 1, j - 1].h_outflow)) + Modelica.Constants.eps*T[i, j, K + 1])/(alpha*abs(heatPort_cei[i - 1, j - 1].Q_flow) + beta*max(0, fluidPort_cei[i - 1,
        j - 1].m_flow) + Modelica.Constants.eps);                                                                                                    //Tmix_floor[i-1,j-1];

      // heat flows
      heatPort_floor[i-1,j-1].Q_flow = Db_M[i-1,j-1,1]*(heatPort_floor[i-1,j-1].T - T[i,j,2]);
      heatPort_cei[i-1,j-1].Q_flow   = Dt_M[i-1,j-1,K]*(heatPort_cei[i-1,j-1].T   - T[i,j,K+1]);

      end for;
  end for;

  annotation (Diagram(graphics), Icon(graphics));
end Room3D_CFDwall_boundaries_old;
