within CFDModelica.Flow3D;
partial model Room3D_CFDwall_walls
  //extends OOMS_BEPI.SubZonalModels_Extension.Flow3D_a.Room3D_CFDwall_boundaries;
  extends Room3D_CFDwall_base;
equation

  //////////////////////////////////////////////////
  // AVERAGE H_coeff
  h_avg_left + T_h_mean*der(h_avg_left) = sum(sum(dy[j]*dz[k]
    *h_left[j, k]/(dx[1]*0.5) for k in 1:
    K) for j in 1:J)/(roomHeight*roomWidth);
  h_avg_right + T_h_mean*der(h_avg_right) = sum(sum(dy[j]*dz[k]
    *h_right[j, k]/(dx[I]*0.5) for k in 1:
    K) for j in 1:J)/(roomHeight*roomWidth);
  h_avg_front + T_h_mean*der(h_avg_front) = sum(sum(dx[i]*dz[k]
    *h_front[i, k]/(dy[1]*0.5) for k in 1:
    K) for i in 1:I)/(roomHeight*roomBase);
  h_avg_rear + T_h_mean*der(h_avg_rear) = sum(sum(dx[i]*dz[k]
    *h_rear[i, k]/(dy[J]*0.5) for k in 1:
    K) for i in 1:I)/(roomHeight*roomBase);
  h_avg_floor + T_h_mean*der(h_avg_floor)= sum( sum(dy[j]*dx[i]*h_floor[i,j]/(dz[1]*0.5) for i in 1:I) for j in 1:J)/(roomBase*roomWidth);
  h_avg_cei   + T_h_mean*der(h_avg_cei)  = sum( sum(dy[j]*dx[i]*h_cei[i,j]/(dz[K]*0.5) for i in 1:I) for j in 1:J)/(roomBase*roomWidth);

  /////////////////////////////////////////////////
  // COMPUTATION OF THE COEFFICIENTS

  // (1,1,1) CORNER
  De_M[1,1,1] = dy[1]*dz[1]*gamma*(1+(MUx[2,2,2] - mu)/mu)/(dx[1]*0.5+dx[2]*0.5);
  Dw_M[1,1,1] = dy[1]*dz[1]/(dx[1]*0.5)*(if Hmean then h_avg_left else h_left[1,1]);
  Dn_M[1,1,1] = dx[1]*dz[1]*gamma*(1+(MUy[2,2,2] - mu)/mu)/(dy[1]*0.5+dy[2]*0.5);
  Ds_M[1,1,1] = dy[1]*dz[1]/(dy[1]*0.5)*(if Hmean then h_avg_front else h_front[1,1]);
  Dt_M[1,1,1] = dx[1]*dy[1]*gamma*(1+(MUz[2,2,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
  Db_M[1,1,1] = dx[1]*dy[1]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[1,1]);
  h_left[1, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[2, 2, 2] + MUz[2, 2, 2])/2,
    dx[1]*0.5,
    Hleft,
    Hfixed);
  h_front[1, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, 2, 2] + MUz[2, 2, 2])/2,
    dy[1]*0.5,
    Hfront,
    Hfixed);
  h_floor[1, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, 2, 2] + MUy[2, 2, 2])/2,
    dz[1]*0.5,
    Hfloor,
    Hfixed);

  // (I,1,1) CORNER
  De_M[I,1,1] = dy[1]*dz[1]/(dx[I]*0.5)*(if Hmean then h_avg_right else h_right[1,1]);
  Dw_M[I,1,1] = dy[1]*dz[1]*gamma*(1+(MUx[I,2,2] - mu)/mu)/(dx[I-1]*0.5+dx[I]*0.5);
  Dn_M[I,1,1] = dx[I]*dz[1]*gamma*(1+(MUy[I+1,2,2] - mu)/mu)/(dy[1]*0.5+dy[2]*0.5);
  Ds_M[I,1,1] = dx[I]*dz[1]/(dy[1]*0.5)*(if Hmean then h_avg_front else h_front[I,1]);
  Dt_M[I,1,1] = dx[I]*dy[1]*gamma*(1+(MUz[I+1,2,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
  Db_M[I,1,1] = dx[I]*dy[1]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[I,1]);
  h_right[1, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[I + 1, 2, 2] + MUz[I + 1, 2, 2])/2,
    dx[I]*0.5,
    Hright,
    Hfixed);
  h_front[I, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, 2, 2] + MUz[I + 1, 2, 2])/2,
    dy[1]*0.5,
    Hfront,
    Hfixed);
  h_floor[I, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, 2, 2] + MUy[I + 1, 2, 2])/2,
    dz[1]*0.5,
    Hfloor,
    Hfixed);

  // (1,1,K) CORNER
  De_M[1,1,K] = dy[1]*dz[K]*gamma*(1+(MUx[2,2,K+1] - mu)/mu)/(dx[1]*0.5+dx[2]*0.5);
  Dw_M[1,1,K] = dy[1]*dz[K]/(dx[1]*0.5)*(if Hmean then h_avg_left else h_left[1,K]);
  Dn_M[1,1,K] = dx[1]*dz[K]*gamma*(1+(MUy[2,2,K+1] - mu)/mu)/(dy[1]*0.5+dy[2]*0.5);
  Ds_M[1,1,K] = dx[1]*dz[K]/(dy[1]*0.5)*(if Hmean then h_avg_front else h_front[1,K]);
  Dt_M[1,1,K] = dx[1]*dy[1]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[1,1]);
  Db_M[1,1,K] = dx[1]*dy[1]*gamma*(1+(MUz[2,2,K] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
  h_left[1, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[2, 2, K + 1] + MUz[2, 2, K])/2,
    dx[1]*0.5,
    Hleft,
    Hfixed);
  h_front[1, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUz[2, 2, K] + MUx[2, 2, K + 1])/2,
    dy[1]*0.5,
    Hfront,
    Hfixed);
  h_cei[1, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, 2, K + 1] + MUy[2, 2, K + 1])/2,
    dz[K]*0.5,
    Hcei,
    Hfixed);

  // (I,1,K) CORNER
  De_M[I,1,K] = dy[1]*dz[K]/(dx[I]*0.5)*(if Hmean then h_avg_right else h_right[1,K]);
  Dw_M[I,1,K] = dy[1]*dz[K]*gamma*(1+(MUx[I,2,K+1] - mu)/mu)/(dx[I-1]*0.5+dx[I]*0.5);
  Dn_M[I,1,K] = dx[I]*dz[K]*gamma*(1+(MUy[I+1,2,K+1] - mu)/mu)/(dy[1]*0.5+dy[2]*0.5);
  Ds_M[I,1,K] = dx[I]*dz[K]/(dy[1]*0.5)*(if Hmean then h_avg_front else h_front[I,K]);
  Dt_M[I,1,K] = dx[I]*dy[1]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[I,1]);
  Db_M[I,1,K] = dx[I]*dy[1]*gamma*(1+(MUz[I+1,2,K] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
  h_right[1, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[I + 1, 2, K + 1] + MUz[I + 1, 2, K])/2,
    dx[I]*0.5,
    Hright,
    Hfixed);
  h_front[I, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, 2, K + 1] + MUz[I + 1, 2, K])/2,
    dy[1]*0.5,
    Hfront,
    Hfixed);
  h_cei[I, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, 2, K + 1] + MUy[I + 1, 2, K + 1])/2,
    dz[K]*0.5,
    Hcei,
    Hfixed);

  // (1,J,1) CORNER
  De_M[1,J,1] = dy[J]*dz[1]*gamma*(1+(MUx[2,J+1,2] - mu)/mu)/(dx[1]*0.5+dx[2]*0.5);
  Dw_M[1,J,1] = dy[J]*dz[1]/(dx[1]*0.5)*(if Hmean then h_avg_left else h_left[J,1]);
  Dn_M[1,J,1] = dx[1]*dz[1]/(dy[J]*0.5)*(if Hmean then h_avg_rear else h_rear[1,1]);
  Ds_M[1,J,1] = dx[1]*dz[1]*gamma*(1+(MUy[2,J,2] - mu)/mu)/(dy[J]*0.5+dy[J-1]*0.5);
  Dt_M[1,J,1] = dx[1]*dy[J]*gamma*(1+(MUz[2,J+1,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
  Db_M[1,J,1] = dx[1]*dy[J]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[1,J]);
  h_left[J, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[2, J, 2] + MUz[2, J + 1, 2])/2,
    dx[1]*0.5,
    Hleft,
    Hfixed);
  h_rear[1, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, J + 1, 2] + MUz[2, J + 1, 2])/2,
    dy[J]*0.5,
    Hrear,
    Hfixed);
  h_floor[1, J] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, J + 1, 2] + MUy[2, J, 2])/2,
    dz[1]*0.5,
    Hfloor,
    Hfixed);

  // (I,J,1) CORNER
  De_M[I,J,1] = dy[J]*dz[1]/(dx[I]*0.5)*(if Hmean then h_avg_right else h_right[J,1]);
  Dw_M[I,J,1] = dy[J]*dz[1]*gamma*(1+(MUx[I,J+1,2] - mu)/mu)/(dx[I-1]*0.5+dx[I]*0.5);
  Dn_M[I,J,1] = dx[I]*dz[1]/(dy[J]*0.5)*(if Hmean then h_avg_rear else h_rear[I,1]);
  Ds_M[I,J,1] = dx[I]*dz[1]*gamma*(1+(MUy[I+1,J,2] - mu)/mu)/(dy[J]*0.5+dy[J-1]*0.5);
  Dt_M[I,J,1] = dx[I]*dy[J]*gamma*(1+(MUz[I+1,J+1,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
  Db_M[I,J,1] = dx[I]*dy[J]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[I,J]);
  h_right[J, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[I + 1, J, 2] + MUz[I + 1, J + 1, 2])/2,
    dx[I]*0.5,
    Hright,
    Hfixed);
  h_rear[I, 1] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, J + 1, 2] + MUz[I + 1, J + 1, 2])/2,
    dy[J]*0.5,
    Hrear,
    Hfixed);
  h_floor[I, J] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, J + 1, 2] + MUy[I + 1, J, 2])/2,
    dz[1]*0.5,
    Hfloor,
    Hfixed);

  // (1,J,K) CORNER
  De_M[1,J,K] = dy[J]*dz[K]*gamma*(1+(MUx[2,J+1,K+1] - mu)/mu)/(dx[1]*0.5+dx[2]*0.5);
  Dw_M[1,J,K] = dy[J]*dz[K]/(dx[1]*0.5)*(if Hmean then h_avg_left else h_left[J,K]);
  Dn_M[1,J,K] = dx[1]*dz[K]/(dy[J]*0.5)*(if Hmean then h_avg_rear else h_rear[1,K]);
  Ds_M[1,J,K] = dx[1]*dz[K]*gamma*(1+(MUy[2,J,K+1] - mu)/mu)/(dy[J]*0.5+dy[J-1]*0.5);
  Dt_M[1,J,K] = dx[1]*dy[J]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[1,J]);
  Db_M[1,J,K] = dx[1]*dy[J]*gamma*(1+(MUz[2,J+1,K] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);
  h_left[J, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[2, J, K + 1] + MUz[2, J + 1, K])/2,
    dx[1]*0.5,
    Hleft,
    Hfixed);
  h_rear[1, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, J + 1, K + 1] + MUz[2, J + 1, K])/2,
    dy[J]*0.5,
    Hrear,
    Hfixed);
  h_cei[1, J] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[2, J + 1, K + 1] + MUy[2, J, K + 1])/2,
    dz[K]*0.5,
    Hcei,
    Hfixed);

  // (I,J,K) CORNER
  De_M[I,J,K] = dy[J]*dz[K]/(dx[I]*0.5)*(if Hmean then h_avg_right else h_right[J,K]);
  Dw_M[I,J,K] = dy[J]*dz[K]*gamma*(1+(MUx[I,J+1,K+1] - mu)/mu)/(dx[I-1]*0.5+dx[I]*0.5);
  Dn_M[I,J,K] = dx[I]*dz[K]/(dy[J]*0.5)*(if Hmean then h_avg_rear else h_rear[I,K]);
  Ds_M[I,J,K] = dx[I]*dz[K]*gamma*(1+(MUy[I+1,J,K+1] - mu)/mu)/(dy[J]*0.5+dy[J-1]*0.5);
  Dt_M[I,J,K] = dx[I]*dy[J]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[I,J]);
  Db_M[I,J,K] = dx[I]*dy[J]*gamma*(1+(MUz[I+1,J+1,K] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);
  h_right[J, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUy[I + 1, J, K + 1] + MUz[I + 1, J + 1, K + 1])/2,
    dx[I]*0.5,
    Hright,
    Hfixed);
  h_rear[I, K] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, J + 1, K + 1] + MUz[I + 1, J + 1, K + 1])/2,
    dy[J]*0.5,
    Hrear,
    Hfixed);
  h_cei[I, J] = CFDModelica.Functions.H_walls(
    gamma,
    mu,
    (MUx[I, J + 1, K + 1] + MUy[I + 1, J, K + 1])/2,
    dz[K]*0.5,
    Hcei,
    Hfixed);

  // X-ORIENTED EDGES
  for i in 2:I-1 loop
    // ROW (:,1,1)
    De_M[i,1,1] = dy[1]*dz[1]*gamma*(1+(MUx[i+1,2,2] - mu)/mu)/(dx[i]*0.5+dx[i+1]*0.5);
    Dw_M[i,1,1] = dy[1]*dz[1]*gamma*(1+(MUx[i,2,2] - mu)/mu)/(dx[i]*0.5+dx[i-1]*0.5);
    Dn_M[i,1,1] = dx[i]*dz[1]*gamma*(1+(MUy[i+1,2,2] - mu)/mu)/(dy[1]*0.5+dy[2]*0.5);
    Ds_M[i,1,1] = dx[i]*dz[1]/(dy[1]*0.5)*(if Hmean then h_avg_front else h_front[i,1]);
    Dt_M[i,1,1] = dx[i]*dy[1]*gamma*(1+(MUz[i+1,2,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
    Db_M[i,1,1] = dx[i]*dy[1]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[i,1]);
    h_front[i, 1] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, 2, 2] + MUx[i, 2, 2] + MUz[i + 1, 2, 2])/3,
      dy[1]*0.5,
      Hfront,
      Hfixed);
    h_floor[i, 1] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, 2, 2] + MUx[i, 2, 2] + MUy[i + 1, 2, 2])/3,
      dz[1]*0.5,
      Hfloor,
      Hfixed);
    // ROW (:,1,K)
    De_M[i,1,K] = dy[1]*dz[K]*gamma*(1+(MUx[i+1,2,K+1] - mu)/mu)/(dx[i]*0.5+dx[i+1]*0.5);
    Dw_M[i,1,K] = dy[1]*dz[K]*gamma*(1+(MUx[i,2,K+1] - mu)/mu)/(dx[i]*0.5+dx[i-1]*0.5);
    Dn_M[i,1,K] = dx[i]*dz[K]*gamma*(1+(MUy[i+1,2,K+1] - mu)/mu)/(dy[1]*0.5+dy[2]*0.5);
    Ds_M[i,1,K] = dx[i]*dz[K]/(dy[1]*0.5)*(if Hmean then h_avg_front else h_front[i,K]);
    Dt_M[i,1,K] = dx[i]*dy[1]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[i,1]);
    Db_M[i,1,K] = dx[i]*dy[1]*gamma*(1+(MUz[i+1,2,K] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);
    h_front[i, K] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, 2, K + 1] + MUx[i, 2, K + 1] + MUz[i + 1, 2, K])/3,
      dy[1]*0.5,
      Hfront,
      Hfixed);
    h_cei[i, 1] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, 2, K + 1] + MUx[i, 2, K + 1] + MUy[i + 1, 2, K + 1])/3,
      dz[K]*0.5,
      Hcei,
      Hfixed);

    // ROW (:,J,1)
    De_M[i,J,1] = dy[J]*dz[1]*gamma*(1+(MUx[i+1,J+1,2] - mu)/mu)/(dx[i]*0.5+dx[i+1]*0.5);
    Dw_M[i,J,1] = dy[J]*dz[1]*gamma*(1+(MUx[i,J+1,2] - mu)/mu)/(dx[i]*0.5+dx[i-1]*0.5);
    Dn_M[i,J,1] = dx[i]*dz[1]/(dy[J]*0.5)*(if Hmean then h_avg_rear else h_rear[i,1]);
    Ds_M[i,J,1] = dx[i]*dz[1]*gamma*(1+(MUy[i+1,J,2] - mu)/mu)/(dy[J]*0.5+dy[J-1]*0.5);
    Dt_M[i,J,1] = dx[i]*dy[J]*gamma*(1+(MUz[i+1,J+1,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
    Db_M[i,J,1] = dx[i]*dy[J]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[i,J]);
    h_rear[i, 1] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, J + 1, 2] + MUx[i, J + 1, 2] + MUz[i + 1, J + 1, 2])/3,
      dy[J]*0.5,
      Hrear,
      Hfixed);
    h_floor[i, J] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, J + 1, 2] + MUx[i, J + 1, 2] + MUy[i + 1, J, 2])/3,
      dz[1]*0.5,
      Hfloor,
      Hfixed);

    // ROW (:,J,K)
    De_M[i,J,K] = dy[J]*dz[K]*gamma*(1+(MUx[i+1,J+1,K+1] - mu)/mu)/(dx[i]*0.5+dx[i+1]*0.5);
    Dw_M[i,J,K] = dy[J]*dz[K]*gamma*(1+(MUx[i,J+1,K+1] - mu)/mu)/(dx[i]*0.5+dx[i-1]*0.5);
    Dn_M[i,J,K] = dx[i]*dz[K]/(dy[J]*0.5)*(if Hmean then h_avg_rear else h_rear[i,K]);
    Ds_M[i,J,K] = dx[i]*dz[K]*gamma*(1+(MUy[i+1,J,K+1] - mu)/mu)/(dy[J]*0.5+dy[J-1]*0.5);
    Dt_M[i,J,K] = dx[i]*dy[J]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[i,J]);
    Db_M[i,J,K] = dx[i]*dy[J]*gamma*(1+(MUz[i+1,J+1,K] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);
    h_rear[i, K] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, J + 1, K + 1] + MUx[i, J + 1, K + 1] + MUz[i + 1, J + 1, K])/
        3,
      dy[J]*0.5,
      Hrear,
      Hfixed);
    h_cei[i, J] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[i + 1, J + 1, K + 1] + MUx[i, J + 1, K + 1] + MUy[i + 1, J, K + 1])/
        3,
      dz[K]*0.5,
      Hcei,
      Hfixed);

  end for;

  // Y-ORIENTED EDGES
  for j in 2:J-1 loop
    // ROW (1,:,1)
    De_M[1,j,1] = dy[j]*dz[1]*gamma*(1+(MUx[2,j+1,2] - mu)/mu)/(dx[1]*0.5 + dx[2]*0.5);
    Dw_M[1,j,1] = dy[j]*dz[1]/(dx[1]*0.5)*(if Hmean then h_avg_left else h_left[j,1]);
    Dn_M[1,j,1] = dx[1]*dz[1]*gamma*(1+(MUy[2,j+1,2] - mu)/mu)/(dy[j]*0.5 + dy[j+1]*0.5);
    Ds_M[1,j,1] = dx[1]*dz[1]*gamma*(1+(MUy[2,j,2] - mu)/mu)/(dy[j]*0.5 + dy[j-1]*0.5);
    Dt_M[1,j,1] = dx[1]*dy[j]*gamma*(1+(MUz[2,j+1,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
    Db_M[1,j,1] = dx[1]*dy[j]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[1,j]);
    h_left[j, 1] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUz[2, j + 1, 2] + MUy[2, j + 1, 2] + MUy[2, j, 2])/3,
      dx[1]*0.5,
      Hleft,
      Hfixed);
    h_floor[1, j] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[2, j + 1, 2] + MUy[2, j + 1, 2] + MUy[2, j, 2])/3,
      dz[1]*0.5,
      Hfloor,
      Hfixed);

    // ROW (1,:,K)
    De_M[1,j,K] = dy[j]*dz[K]*gamma*(1+(MUx[2,j+1,K+1] - mu)/mu)/(dx[1]*0.5 + dx[2]*0.5);
    Dw_M[1,j,K] = dy[j]*dz[K]/(dx[1]*0.5)*(if Hmean then h_avg_left else h_left[j,K]);
    Dn_M[1,j,K] = dx[1]*dz[K]*gamma*(1+(MUy[2,j+1,K+1] - mu)/mu)/(dy[j]*0.5 + dy[j+1]*0.5);
    Ds_M[1,j,K] = dx[1]*dz[K]*gamma*(1+(MUy[2,j,K+1] - mu)/mu)/(dy[j]*0.5 + dy[j-1]*0.5);
    Dt_M[1,j,K] = dx[1]*dy[j]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[1,j]);
    Db_M[1,j,K] = dx[1]*dy[j]*gamma*(1+(MUz[2,j+1,K+1] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);
    h_left[j, K] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUz[2, j + 1, K + 1] + MUy[2, j + 1, K + 1] + MUy[2, j, K + 1])/3,
      dx[1]*0.5,
      Hleft,
      Hfixed);
    h_cei[1, j] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[2, j + 1, K + 1] + MUy[2, j + 1, K + 1] + MUy[2, j, K + 1])/3,
      dz[K]*0.5,
      Hcei,
      Hfixed);

    // ROW (I,:,1)
    De_M[I,j,1] = dy[j]*dz[1]/(dx[I]*0.5)*(if Hmean then h_avg_right else h_right[j,1]);
    Dw_M[I,j,1] = dy[j]*dz[1]*gamma*(1+(MUx[I,j+1,2] - mu)/mu)/(dx[I]*0.5 + dx[I-1]*0.5);
    Dn_M[I,j,1] = dx[I]*dz[1]*gamma*(1+(MUy[I+1,j+1,K+1] - mu)/mu)/(dy[j]*0.5 + dy[j+1]*0.5);
    Ds_M[I,j,1] = dx[I]*dz[1]*gamma*(1+(MUy[I+1,j,K+1] - mu)/mu)/(dy[j]*0.5 + dy[j-1]*0.5);
    Dt_M[I,j,1] = dx[I]*dy[j]*gamma*(1+(MUz[I+1,j+1,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
    Db_M[I,j,1] = dx[I]*dy[j]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[I,j]);
    h_right[j, 1] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUy[I + 1, j + 1, K + 1] + MUy[I + 1, j, K + 1] + MUz[I + 1, j + 1, 2])/
        3,
      dx[I]*0.5,
      Hright,
      Hfixed);
    h_floor[I, j] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[I, j + 1, 2] + MUy[I + 1, j + 1, K + 1] + MUy[I + 1, j, K + 1])/3,
      dz[1]*0.5,
      Hfloor,
      Hfixed);

    // ROW (I,:,K)
    De_M[I,j,K] = dy[j]*dz[K]/(dx[I]*0.5)*(if Hmean then h_avg_right else h_right[j,K]);
    Dw_M[I,j,K] = dy[j]*dz[K]*gamma*(1+(MUx[I,j+1,K+1] - mu)/mu)/(dx[I]*0.5 + dx[I-1]*0.5);
    Dn_M[I,j,K] = dx[I]*dz[1]*gamma*(1+(MUy[I+1,j+1,K+1] - mu)/mu)/(dy[j]*0.5 + dy[j+1]*0.5);
    Ds_M[I,j,K] = dx[I]*dz[1]*gamma*(1+(MUy[I+1,j,K+1] - mu)/mu)/(dy[j]*0.5 + dy[j-1]*0.5);
    Dt_M[I,j,K] = dx[I]*dy[j]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[I,j]);
    Db_M[I,j,K] = dx[I]*dy[j]*gamma*(1+(MUz[I+1,j+1,K+1] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);
    h_right[j, K] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUy[I + 1, j + 1, K + 1] + MUy[I + 1, j, K + 1] + MUz[I + 1, j + 1, K +
        1])/3,
      dx[I]*0.5,
      Hright,
      Hfixed);
    h_cei[I, j] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[I, j + 1, K + 1] + MUy[I + 1, j + 1, K + 1] + MUy[I + 1, j, K + 1])/
        3,
      dz[K]*0.5,
      Hcei,
      Hfixed);

  end for;

  // Z-ORIENTED EDGES
  for k in 2:K - 1 loop
    // ROW (1,1,:)
    De_M[1, 1, k] = dy[1]*dz[k]*gamma*(1 +
      (MUx[2, 2, k + 1] - mu)/mu)/(dx[1]*0.5 + dx[2]*0.5);
    Dw_M[1, 1, k] = dy[1]*dz[k]/(dx[1]*0.5)
      *(if Hmean then h_avg_left else h_left[1, k]);
    Dn_M[1, 1, k] = dx[1]*dz[k]*gamma*(1 +
      (MUy[2, 2, k + 1] - mu)/mu)/(dy[1]*0.5 + dy[2]*0.5);
    Ds_M[1, 1, k] = dx[1]*dz[k]/(dy[1]*0.5)
      *(if Hmean then h_avg_front else h_front[1, k]);
    Dt_M[1, 1, k] = dx[1]*dy[1]*gamma*(1 + (MUz[2, 2,
      k + 1] - mu)/mu)/(dz[k]*0.5 + dz[
      k + 1]*0.5);
    Db_M[1, 1, k] = dx[1]*dy[1]*gamma*(1 + (MUz[2, 2,
      k] - mu)/mu)/(dz[k]*0.5 + dz[
      k - 1]*0.5);
    h_left[1, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUy[2, 2, k + 1] + MUz[2, 2, k + 1] + MUz[2, 2, k])/3,
      dx[1]*0.5,
      Hleft,
      Hfixed);
    h_front[1, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[2, 2, k + 1] + MUz[2, 2, k + 1] + MUz[2, 2, k])/3,
      dy[1]*0.5,
      Hfront,
      Hfixed);

    // ROW (1,J,:)
    De_M[1, J, k] = dy[J]*dz[k]*gamma*(1 +
      (MUx[2, J + 1, k + 1] - mu)/mu)/(dx[1]*0.5 + dx[2]*0.5);
    Dw_M[1, J, k] = dy[J]*dz[k]/(dx[1]*0.5)
      *(if Hmean then h_avg_left else h_left[J, k]);
    Dn_M[1, J, k] = dx[1]*dz[k]/(dy[J]*0.5)
      *(if Hmean then h_avg_rear else h_rear[1, k]);
    Ds_M[1, J, k] = dx[1]*dz[k]*gamma*(1 +
      (MUy[2, J, k + 1] - mu)/mu)/(dy[J]*0.5 + dy[J - 1]*0.5);
    Dt_M[1, J, k] = dx[1]*dy[J]*gamma*(1 + (MUz[2, J + 1,
      k + 1] - mu)/mu)/(dz[k]*0.5 + dz[
      k + 1]*0.5);
    Db_M[1, J, k] = dx[1]*dy[J]*gamma*(1 + (MUz[2, J + 1,
      k] - mu)/mu)/(dz[k]*0.5 + dz[
      k - 1]*0.5);
    h_left[J, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUy[2, J, k + 1] + MUz[2, J + 1, k + 1] + MUz[2, J + 1, k])/3,
      dx[1]*0.5,
      Hleft,
      Hfixed);
    h_rear[1, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[2, J + 1, k + 1] + MUz[2, J + 1, k + 1] + MUz[2, J + 1, k])/3,
      dy[J]*0.5,
      Hrear,
      Hfixed);

    // ROW (I,1,:)
    De_M[I, 1, k] = dy[1]*dz[k]/(dx[I]*0.5)
      *(if Hmean then h_avg_right else h_right[1, k]);
    Dw_M[I, 1, k] = dy[1]*dz[k]*gamma*(1 +
      (MUx[I, 2, k + 1] - mu)/mu)/(dx[I]*0.5 + dx[I - 1]*0.5);
    Dn_M[I, 1, k] = dx[I]*dz[k]*gamma*(1 +
      (MUy[I + 1, 2, k + 1] - mu)/mu)/(dy[1]*0.5 + dy[2]*0.5);
    Ds_M[I, 1, k] = dx[I]*dz[k]/(dy[1]*0.5)
      *(if Hmean then h_avg_front else h_front[I, k]);
    Dt_M[I, 1, k] = dx[I]*dy[1]*gamma*(1 + (MUz[I + 1, 2,
      k + 1] - mu)/mu)/(dz[k]*0.5 + dz[
      k + 1]*0.5);
    Db_M[I, 1, k] = dx[I]*dy[1]*gamma*(1 + (MUz[I + 1, 2,
      k] - mu)/mu)/(dz[k]*0.5 + dz[
      k - 1]*0.5);
    h_right[1, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUy[I + 1, 2, k + 1] + MUz[I + 1, 2, k + 1] + MUz[I + 1, 2, k])/3,
      dx[I]*0.5,
      Hright,
      Hfixed);
    h_front[I, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[I, 2, k + 1] + MUz[I + 1, 2, k + 1] + MUz[I + 1, 2, k])/3,
      dy[1]*0.5,
      Hfront,
      Hfixed);

    // ROW (I,J,:)
    De_M[I, J, k] = dy[1]*dz[k]/(dx[I]*0.5)
      *(if Hmean then h_avg_right else h_right[J, k]);
    Dw_M[I, J, k] = dy[1]*dz[k]*gamma*(1 +
      (MUx[I, J + 1, k + 1] - mu)/mu)/(dx[I]*0.5 + dx[I - 1]*
      0.5);
    Dn_M[I, J, k] = dx[I]*dz[k]/(dy[J]*0.5)
      *(if Hmean then h_avg_rear else h_rear[I, k]);
    Ds_M[I, J, k] = dx[I]*dz[k]*gamma*(1 +
      (MUy[I + 1, J, k + 1] - mu)/mu)/(dy[J]*0.5 + dy[J - 1]*
      0.5);
    Dt_M[I, J, k] = dx[I]*dy[J]*gamma*(1 + (MUz[I + 1, J + 1,
      k + 1] - mu)/mu)/(dz[k]*0.5 + dz[
      k + 1]*0.5);
    Db_M[I, J, k] = dx[I]*dy[J]*gamma*(1 + (MUz[I + 1, J + 1,
      k] - mu)/mu)/(dz[k]*0.5 + dz[
      k - 1]*0.5);
    h_right[J, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUy[I + 1, J, k + 1] + MUz[I + 1, J + 1, k + 1] + MUz[I + 1, J + 1, k])/
        3,
      dx[I]*0.5,
      Hright,
      Hfixed);
    h_rear[I, k] = CFDModelica.Functions.H_walls(
      gamma,
      mu,
      (MUx[I, J + 1, k + 1] + MUz[I + 1, J + 1, k + 1] + MUz[I + 1, J + 1, k])/
        3,
      dy[J]*0.5,
      Hrear,
      Hfixed);

  end for;

  // RIGHT AND LEFT WALLS
  for j in 2:J-1 loop
    for k in 2:K - 1 loop
      // right
      De_M[I, j, k] = dy[j]*dz[k]/(dx[I]*0.5)
        *(if Hmean then h_avg_right else h_right[j, k]);
      Dw_M[I, j, k] = dy[j]*dz[k]*gamma*(1
         + (MUx[I, j + 1, k + 1] - mu)/mu)/(dx[I]*0.5 + dx[I -
        1]*0.5);
      Dn_M[I, j, k] = dz[k]*dx[I]*gamma*(1
         + (MUy[I + 1, j + 1, k + 1] - mu)/mu)/(dy[j]*0.5 +
        dy[j + 1]*0.5);
      Ds_M[I, j, k] = dz[k]*dx[I]*gamma*(1
         + (MUy[I + 1, j, k + 1] - mu)/mu)/(dy[j]*0.5 + dy[j -
        1]*0.5);
      Dt_M[I, j, k] = dy[j]*dx[I]*gamma*(1 + (MUz[I + 1, j +
        1, k + 1] - mu)/mu)/(dz[k]*0.5 +
        dz[k + 1]*0.5);
      Db_M[I, j, k] = dy[j]*dx[I]*gamma*(1 + (MUz[I + 1, j +
        1, k] - mu)/mu)/(dz[k]*0.5 + dz[
        k - 1]*0.5);

      h_right[j, k] = CFDModelica.Functions.H_walls(
        gamma,
        mu,
        (MUz[I + 1, j + 1, k] + MUz[I + 1, j + 1, k + 1] + MUy[I + 1, j, k + 1]
           + MUy[I + 1, j + 1, k + 1])/4,
        dx[I]*0.5,
        Hright,
        Hfixed);

      // left
      De_M[1, j, k] = dy[j]*dz[k]*gamma*(1
         + (MUx[2, j + 1, k + 1] - mu)/mu)/(dx[1]*0.5 + dx[2]
        *0.5);
      Dw_M[1, j, k] = dy[j]*dz[k]/(dx[1]*0.5)
        *(if Hmean then h_avg_left else h_left[j, k]);
      Dn_M[1, j, k] = dz[k]*dx[1]*gamma*(1
         + (MUy[2, j + 1, k + 1] - mu)/mu)/(dy[j]*0.5 + dy[j +
        1]*0.5);
      Ds_M[1, j, k] = dz[k]*dx[1]*gamma*(1
         + (MUy[2, j, k + 1] - mu)/mu)/(dy[j]*0.5 + dy[j - 1]
        *0.5);
      Dt_M[1, j, k] = dy[j]*dx[1]*gamma*(1 + (MUz[2, j + 1,
        k + 1] - mu)/mu)/(dz[k]*0.5 + dz[
        k + 1]*0.5);
      Db_M[1, j, k] = dy[j]*dx[1]*gamma*(1 + (MUz[2, j + 1,
        k] - mu)/mu)/(dz[k]*0.5 + dz[
        k - 1]*0.5);

      h_left[j, k] = CFDModelica.Functions.H_walls(
        gamma,
        mu,
        (MUz[2, j + 1, k] + MUz[2, j + 1, k + 1] + MUy[2, j, k + 1] + MUy[2, j
           + 1, k + 1])/4,
        dx[1]*0.5,
        Hleft,
        Hfixed);
    end for;
  end for;
  // FRONT AND REAR WALLS
  for i in 2:I-1 loop
    for k in 2:K - 1 loop
      // front
      De_M[i, 1, k] = dy[1]*dz[k]*gamma*(1
         + (MUx[i + 1, 2, k + 1] - mu)/mu)/(dx[i]*0.5 + dx[i +
        1]*0.5);
      Dw_M[i, 1, k] = dy[1]*dz[k]*gamma*(1
         + (MUx[i, 2, k + 1] - mu)/mu)/(dx[i]*0.5 + dx[i - 1]
        *0.5);
      Dn_M[i, 1, k] = dz[k]*dx[i]*gamma*(1
         + (MUy[i + 1, 2, k + 1] - mu)/mu)/(dy[1]*0.5 + dy[2]
        *0.5);
      Ds_M[i, 1, k] = dz[k]*dx[i]/(dy[1]*0.5)
        *(if Hmean then h_avg_front else h_front[i, k]);
      Dt_M[i, 1, k] = dy[1]*dx[i]*gamma*(1 + (MUz[i + 1, 2,
        k + 1] - mu)/mu)/(dz[k]*0.5 + dz[
        k + 1]*0.5);
      Db_M[i, 1, k] = dy[1]*dx[i]*gamma*(1 + (MUz[i + 1, 2,
        k] - mu)/mu)/(dz[k]*0.5 + dz[
        k - 1]*0.5);

      h_front[i, k] = CFDModelica.Functions.H_walls(
        gamma,
        mu,
        (MUz[i + 1, 2, k] + MUz[i + 1, 2, k + 1] + MUx[i + 1, 2, k + 1] + MUx[i,
          2, k + 1])/4,
        dy[1]*0.5,
        Hfront,
        Hfixed);

      // rear
      De_M[i, J, k] = dy[J]*dz[k]*gamma*(1
         + (MUx[i + 1, J + 1, k + 1] - mu)/mu)/(dx[i]*0.5 +
        dx[i + 1]*0.5);
      Dw_M[i, J, k] = dy[J]*dz[k]*gamma*(1
         + (MUx[i, J + 1, k + 1] - mu)/mu)/(dx[i]*0.5 + dx[i -
        1]*0.5);
      Dn_M[i, J, k] = dz[k]*dx[i]/(dy[J]*0.5)
        *(if Hmean then h_avg_rear else h_rear[i, k]);
      Ds_M[i, J, k] = dz[k]*dx[i]*gamma*(1
         + (MUy[i + 1, J + 1, k + 1] - mu)/mu)/(dy[J - 1]*0.5
         + dy[J]*0.5);
      Dt_M[i, J, k] = dy[J]*dx[i]*gamma*(1 + (MUz[i + 1, J +
        1, k + 1] - mu)/mu)/(dz[k]*0.5 +
        dz[k + 1]*0.5);
      Db_M[i, J, k] = dy[J]*dx[i]*gamma*(1 + (MUz[i + 1, J +
        1, k] - mu)/mu)/(dz[k]*0.5 + dz[
        k - 1]*0.5);

      h_rear[i, k] = CFDModelica.Functions.H_walls(
        gamma,
        mu,
        (MUz[i + 1, J + 1, k] + MUz[i + 1, J + 1, k + 1] + MUx[i + 1, J + 1, k
           + 1] + MUx[i, J + 1, k + 1])/4,
        dy[J]*0.5,
        Hrear,
        Hfixed);
    end for;
  end for;
  // CEILING AND FLOOR
  for i in 2:I-1 loop
    for j in 2:J-1 loop
      // ceiling
      De_M[i,j,K] = dy[j]*dz[K]*gamma*(1+(MUx[i+1,j+1,K+1] - mu)/mu)/(dx[i]*0.5+dx[i+1]*0.5);
      Dw_M[i,j,K] = dy[j]*dz[K]*gamma*(1+(MUx[i,j+1,K+1] - mu)/mu)/(dx[i]*0.5+dx[i-1]*0.5);
      Dn_M[i,j,K] = dx[i]*dz[K]*gamma*(1+(MUy[i+1,j+1,K+1] - mu)/mu)/(dy[j]*0.5+dy[j+1]*0.5);
      Ds_M[i,j,K] = dx[i]*dz[K]*gamma*(1+(MUy[i+1,j,K+1] - mu)/mu)/(dy[j]*0.5+dy[j-1]*0.5);
      Dt_M[i,j,K] = dy[j]*dx[i]/(dz[K]*0.5)*(if Hmean then h_avg_cei else h_cei[i,j]);
      Db_M[i,j,K] = dy[j]*dx[i]*gamma*(1+(MUz[i+1,j+1,K] - mu)/mu)/(dz[K]*0.5+dz[K-1]*0.5);

      h_cei[i, j] = CFDModelica.Functions.H_walls(
        gamma,
        mu,
        (MUx[i, j + 1, K + 1] + MUx[i + 1, j + 1, K + 1] + MUy[i + 1, j, K + 1]
           + MUy[i + 1, j + 1, K + 1])/4,
        dz[K]*0.5,
        Hcei,
        Hfixed);
      // floor
      De_M[i,j,1] = dy[j]*dz[1]*gamma*(1+(MUx[i+1,j+1,2] - mu)/mu)/(dx[i]*0.5+dx[i+1]*0.5);
      Dw_M[i,j,1] = dy[j]*dz[1]*gamma*(1+(MUx[i,j+1,2] - mu)/mu)/(dx[i]*0.5+dx[i-1]*0.5);
      Dn_M[i,j,1] = dx[i]*dz[1]*gamma*(1+(MUy[i+1,j+1,2] - mu)/mu)/(dy[j]*0.5+dy[j+1]*0.5);
      Ds_M[i,j,1] = dx[i]*dz[1]*gamma*(1+(MUy[i+1,j,2] - mu)/mu)/(dy[j]*0.5+dy[j-1]*0.5);
      Dt_M[i,j,1] = dy[j]*dx[i]*gamma*(1+(MUz[i+1,j+1,2] - mu)/mu)/(dz[1]*0.5+dz[2]*0.5);
      Db_M[i,j,1] = dy[j]*dx[i]/(dz[1]*0.5)*(if Hmean then h_avg_floor else h_floor[i,j]);

      h_floor[i, j] = CFDModelica.Functions.H_walls(
        gamma,
        mu,
        (MUx[i, j + 1, 2] + MUx[i + 1, j + 1, 2] + MUy[i + 1, j, 2] + MUy[i + 1,
          j + 1, 2])/4,
        dz[1]*0.5,
        Hfloor,
        Hfixed);
    end for;
  end for;

  // INTERIOR
  for i in 2:I-1 loop
    for j in 2:J-1 loop
      for k in 2:K - 1 loop
        De_M[i, j, k] = dy[j]*dz[k]*gamma*
          (1 + (MUx[i + 1, j + 1, k + 1] - mu)/mu)/(dx[i]*0.5
           + dx[i + 1]*0.5);
        Dw_M[i, j, k] = dy[j]*dz[k]*gamma*
          (1 + (MUx[i, j + 1, k + 1] - mu)/mu)/(dx[i]*0.5 +
          dx[i - 1]*0.5);

        Dn_M[i, j, k] = dx[i]*dz[k]*gamma*
          (1 + (MUy[i + 1, j + 1, k + 1] - mu)/mu)/(dy[j]*0.5
           + dy[j + 1]*0.5);
        Ds_M[i, j, k] = dx[i]*dz[k]*gamma*
          (1 + (MUy[i + 1, j, k + 1] - mu)/mu)/(dy[j]*0.5 +
          dy[j - 1]*0.5);

        Dt_M[i, j, k] = dy[j]*dx[i]*gamma*(1 + (MUz[i + 1, j +
          1, k + 1] - mu)/mu)/(dz[k]*0.5 +
          dz[k + 1]*0.5);
        Db_M[i, j, k] = dy[j]*dx[i]*gamma*(1 + (MUz[i + 1, j +
          1, k] - mu)/mu)/(dz[k]*0.5 + dz[
          k - 1]*0.5);
      end for;
    end for;
  end for;

end Room3D_CFDwall_walls;
