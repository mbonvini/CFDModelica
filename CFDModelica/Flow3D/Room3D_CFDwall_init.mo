within CFDModelica.Flow3D;
partial model Room3D_CFDwall_init
  extends Room3D_CFDwall_base;
initial equation

  // initial conditions for pressures and densities inside the volumes
  for i in 2:I+1 loop
    for j in 2:J+1 loop
      for k in 2:K + 1 loop

        /************************************************************************************************
        *********** implicit initalisation **************************************************************
        *************************************************************************************************
        rho[i,j,k] = rho_o + ComprCoeff*(rho[i,j,k]*g*height[i-1,j-1,k-1]) - ThermalExpCoeff*(Tstart-To);
        P[i-1,j-1,k-1] = rho[i,j,k]*g*height[i-1,j-1,k-1] + Medium.p0 ;
        *************************************************************************************************/

        /************************************************************************************************
        *********** implicit initalisation - WITH MEDIA *************************************************
        ************************************************************************************************/
        if energy then
          if not gravity then
            rho[i, j, k] = Medium.density_pT(0, Tstart);
          else
            rho[i, j, k] = Medium.density_pT(rho[i, j,
              k]*g*height[i - 1, j - 1, k -
              1], Tstart);
            P[i - 1, j - 1, k - 1] = rho[i, j, k]
              *g*height[i - 1, j - 1, k - 1];
          end if;
        else
          if not gravity then
            P[i - 1, j - 1, k - 1] = 0;
          else
            rho[i, j, k] = Medium.density_pT(rho[i, j,
              k]*g*height[i - 1, j - 1, k -
              1], Tstart);
            P[i - 1, j - 1, k - 1] = rho[i, j, k]
              *g*height[i - 1, j - 1, k - 1];
          end if;
        end if;

        /********** explicit initialisation *************************************************************/
        //rho[i,j,k] = (ThermalExpCoeff*Tstart - ThermalExpCoeff*To - rho_o)/(ComprCoeff*g*height[i-1,j-1,k-1] - 1);
        //P[i-1,j-1,k-1] = 0*Medium.p0 + (g*height[i-1,j-1,k-1]*ThermalExpCoeff*Tstart-g*height[i-1,j-1,k-1]*ThermalExpCoeff*To-g*height[i-1,j-1,k-1]*rho_o)/(ComprCoeff*g*height[i-1,j-1,k-1]-1);
        /************************************************************************************************/

      end for;
    end for;
  end for;

  // initial value for variable and averaged convective heat transfer coefficient
  h_avg_left = gamma*2/dx[1];
  h_avg_right = gamma*2/dx[I];
  h_avg_front = gamma*2/dy[1];
  h_avg_rear = gamma*2/dy[J];
  h_avg_floor = gamma*2/dz[1];
  h_avg_cei = gamma*2/dz[K];

equation

  // total mass and energy
  Mtot = sum(m);
  Etot = sum(sum(sum(m[i, j, k]*e[i, j, k]
    for i in 1:I) for j in 1:J) for k in 1:K);

  // heat flows through the surfaces
  FluxLeft = sum(sum(Dw_M[1, j, k]*(T[1, j + 1, k
     + 1] - T[2, j + 1, k + 1]) for k in 1
    :K) for j in 1:J);
  FluxRight = sum(sum(De_M[I, j, k]*(T[I + 2, j + 1,
    k + 1] - T[I + 1, j + 1, k + 1]) for
    k in 1:K) for j in 1:J);
  FluxFront = sum(sum(Ds_M[i, 1, k]*(T[i + 1, 1, k
     + 1] - T[i + 1, 2, k + 1]) for k in 1
    :K) for i in 1:I);
  FluxRear = sum(sum(Dn_M[i, J, k]*(T[i + 1, J + 2, k
     + 1] - T[i + 1, J + 1, k + 1]) for k in
        1:K) for i in 1:I);
  FluxFloor   = sum( sum(Db_M[i,j,1]*(T[i+1,j+1,1]  - T[i+1,j+1,2])   for i in 1:I)  for j in 1:J);
  FluxCeiling = sum( sum(Dt_M[i,j,K]*(T[i+1,j+1,K+2]- T[i+1,j+1,K+1]) for i in 1:I)  for j in 1:J);

end Room3D_CFDwall_init;
