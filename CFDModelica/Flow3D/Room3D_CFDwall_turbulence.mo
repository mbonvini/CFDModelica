within CFDModelica.Flow3D;
partial model Room3D_CFDwall_turbulence
  //extends OOMS_BEPI.SubZonalModels_Extension.Flow3D_a.Room3D_CFDwall_massANDenergy;
  extends Room3D_CFDwall_base;
equation

  ////////////////////////////////////////
  // DYAMIC VISCOSITY
  // x-diretion
  ////////////////////////////////////////
  // values on the boundaries
  for i in 1:I+1 loop
    for j in 1:J+2 loop
      // value on the floor
      MUx[i, j, 1] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vx[i, j, 1]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceX_3D(
        I,
        J,
        K,
        i,
        j,
        1,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
      // value on the ceiling
      MUx[i, j, K + 2] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vx[i, j, K + 2]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceX_3D(
        I,
        J,
        K,
        i,
        j,
        K + 2,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
    end for;

    for k in 2:K + 1 loop
      // value on the front
      MUx[i, 1, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vx[i, 1, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceX_3D(
        I,
        J,
        K,
        i,
        1,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
      // value on the rear
      MUx[i, J + 2, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vx[i, J + 2, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceX_3D(
        I,
        J,
        K,
        i,
        J + 2,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
    end for;
  end for;

  // boundary layers
  for i in 1:I+1 loop
    for j in 2:J+1 loop
      // floor
      MUx[i, j, 2] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vx[i, j, 2],
        rho_o,
        CFDModelica.Functions.wallDistanceX_3D(
          I,
          J,
          K,
          i,
          j,
          2,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      // ceiling
      MUx[i, j, K + 1] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vx[i, j, K + 1],
        rho_o,
        CFDModelica.Functions.wallDistanceX_3D(
          I,
          J,
          K,
          i,
          j,
          K + 1,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
    end for;

    for k in 3:K loop
      // front
      MUx[i, 2, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vx[i, 2, k],
        rho_o,
        CFDModelica.Functions.wallDistanceX_3D(
          I,
          J,
          K,
          i,
          2,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      // rear
      MUx[i, J + 1, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vx[i, J + 1, k],
        rho_o,
        CFDModelica.Functions.wallDistanceX_3D(
          I,
          J,
          K,
          i,
          J + 1,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
    end for;
  end for;

  // internal values
  for i in 1:I+1 loop
    for j in 3:J loop
      for k in 3:K loop
        MUx[i, j, k] = if laminar then mu else (mu + Kturb*
          CFDModelica.Functions.sqrtReg(Vx[i, j, k]^2, 1e-8)*rho_o*
          CFDModelica.Functions.wallDistanceX_3D(
          I,
          J,
          K,
          i,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      end for;
    end for;
  end for;

  ////////////////////////////////////////
  // DYAMIC VISCOSITY
  // y-diretion
  ////////////////////////////////////////
  // values on the boundaries
  for j in 1:J+1 loop
    for i in 1:I+2 loop
      // value on the floor
      MUy[i, j, 1] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vy[i, j, 1]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceY_3D(
        I,
        J,
        K,
        i,
        j,
        1,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
      // value on the ceiling
      MUy[i, j, K + 2] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vy[i, j, K + 2]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceY_3D(
        I,
        J,
        K,
        i,
        j,
        K + 2,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
    end for;

    for k in 2:K + 1 loop
      // value on the left
      MUy[1, j, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vy[1, j, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceY_3D(
        I,
        J,
        K,
        1,
        j,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
      // value on the right
      MUy[I + 2, j, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vy[I + 2, j, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceY_3D(
        I,
        J,
        K,
        I + 2,
        j,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
    end for;

  end for;

  // boundary layers
  for j in 1:J+1 loop
    for i in 2:I+1 loop
      // floor
      MUy[i, j, 2] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vy[i, j, 2],
        rho_o,
        CFDModelica.Functions.wallDistanceY_3D(
          I,
          J,
          K,
          i,
          j,
          2,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      // ceiling
      MUy[i, j, K + 1] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vy[i, j, K + 1],
        rho_o,
        CFDModelica.Functions.wallDistanceY_3D(
          I,
          J,
          K,
          i,
          j,
          K + 1,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
    end for;

    for k in 3:K loop
      // left
      MUy[2, j, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vy[2, j, k],
        rho_o,
        CFDModelica.Functions.wallDistanceY_3D(
          I,
          J,
          K,
          2,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      // right
      MUy[I + 1, j, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vy[I + 1, j, k],
        rho_o,
        CFDModelica.Functions.wallDistanceY_3D(
          I,
          J,
          K,
          I + 1,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
    end for;

  end for;

  // internal values
  for j in 1:J+1 loop
    for i in 3:I loop
      for k in 3:K loop
        MUy[i, j, k] = if laminar then mu else (mu + Kturb*
          CFDModelica.Functions.sqrtReg(Vy[i, j, k]^2, 1e-8)*rho_o*
          CFDModelica.Functions.wallDistanceY_3D(
          I,
          J,
          K,
          i,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      end for;
    end for;
  end for;

  ////////////////////////////////////////
  // DYAMIC VISCOSITY
  // z-diretion
  ////////////////////////////////////////
  // values on the boundaries
  for k in 1:K + 1 loop
    for i in 1:I+2 loop
      // value on the front
      MUz[i, 1, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vz[i, 1, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceZ_3D(
        I,
        J,
        K,
        i,
        1,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
      // value on the rear
      MUz[i, J + 2, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vz[i, J + 2, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceZ_3D(
        I,
        J,
        K,
        i,
        J + 2,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
    end for;

    for j in 2:J+1 loop
      // value on the left
      MUz[1, j, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vz[1, j, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceZ_3D(
        I,
        J,
        K,
        1,
        j,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
      // value on the right
      MUz[I + 2, j, k] = if laminar then mu else (mu + Kturb*
        CFDModelica.Functions.sqrtReg(Vz[I + 2, j, k]^2, 1e-8)*rho_o*
        CFDModelica.Functions.wallDistanceZ_3D(
        I,
        J,
        K,
        I + 2,
        j,
        k,
        lengthWall,
        X_frac[:],
        Y_frac[:],
        Z_frac[:],
        roomBase,
        roomWidth,
        roomHeight));
    end for;

  end for;

  // boundary layers
  for k in 1:K + 1 loop
    for i in 2:I+1 loop
      // front
      MUz[i, 2, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vz[i, 2, k],
        rho_o,
        CFDModelica.Functions.wallDistanceZ_3D(
          I,
          J,
          K,
          i,
          2,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      // rear
      MUz[i, J + 1, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vz[i, J + 1, k],
        rho_o,
        CFDModelica.Functions.wallDistanceZ_3D(
          I,
          J,
          K,
          i,
          J + 1,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
    end for;

    for j in 3:J loop
      // left
      MUz[2, j, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vz[2, j, k],
        rho_o,
        CFDModelica.Functions.wallDistanceZ_3D(
          I,
          J,
          K,
          2,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      // right
      MUz[I + 1, j, k] = if laminar then mu else CFDModelica.Functions.MU_walls(
        mu,
        Kturb,
        Vz[I + 1, j, k],
        rho_o,
        CFDModelica.Functions.wallDistanceZ_3D(
          I,
          J,
          K,
          I + 1,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
    end for;

  end for;

  // internal values
  for k in 1:K + 1 loop
    for i in 3:I loop
      for j in 3:J loop
        MUz[i, j, k] = if laminar then mu else (mu + Kturb*
          CFDModelica.Functions.sqrtReg(Vz[i, j, k]^2, 1e-8)*rho_o*
          CFDModelica.Functions.wallDistanceZ_3D(
          I,
          J,
          K,
          i,
          j,
          k,
          lengthWall,
          X_frac[:],
          Y_frac[:],
          Z_frac[:],
          roomBase,
          roomWidth,
          roomHeight));
      end for;
    end for;
  end for;

end Room3D_CFDwall_turbulence;
