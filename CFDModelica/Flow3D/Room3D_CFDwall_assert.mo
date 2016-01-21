within CFDModelica.Flow3D;
partial model Room3D_CFDwall_assert
  extends Room3D_CFDwall_base;
equation
  // control the volum sizes
  assert(abs(sum(X_frac[:])-1)<=1e-6,"error in the x-fraction vector");
  assert(abs(sum(Y_frac[:])-1)<=1e-6,"error in the y-fraction vector");
  assert(abs(sum(Z_frac[:])-1)<=1e-6,"error in the z-fraction vector");

  // control the connectors: thermal and fluid
  // if a thermal connector is connected, the relative thermal one has tobe disconnected
  for j in 1:J loop
    for k in 1:K loop
      assert(abs(fluidPort_left[j, k].m_flow) < Modelica.Constants.eps or abs(
        heatPort_left[j, k].Q_flow) < Modelica.Constants.eps, "(left) there is at least one couple of fluid and thermal connectors coupled");
    end for;
  end for;
  for j in 1:J loop
    for k in 1:K loop
      assert(abs(fluidPort_right[j, k].m_flow) < Modelica.Constants.eps or abs(
        heatPort_right[j, k].Q_flow) < Modelica.Constants.eps, "(right) there is at least one couple of fluid and thermal connectors coupled");
    end for;
  end for;
  for i in 1:I loop
    for k in 1:K loop
      assert(abs(fluidPort_front[i, k].m_flow) < Modelica.Constants.eps or abs(
        heatPort_front[i, k].Q_flow) < Modelica.Constants.eps, "(front) there is at least one couple of fluid and thermal connectors coupled");
    end for;
  end for;
  for i in 1:I loop
    for k in 1:K loop
      assert(abs(fluidPort_rear[i, k].m_flow) < Modelica.Constants.eps or abs(
        heatPort_rear[i, k].Q_flow) < Modelica.Constants.eps, "(rear) there is at least one couple of fluid and thermal connectors coupled");
    end for;
  end for;
  for i in 1:I loop
    for j in 1:J loop
      assert(abs(fluidPort_floor[i, j].m_flow) < Modelica.Constants.eps or abs(
        heatPort_floor[i, j].Q_flow) < Modelica.Constants.eps, "(floor) there is at least one couple of fluid and thermal connectors coupled");
    end for;
  end for;
  for i in 1:I loop
    for j in 1:J loop
      assert(abs(fluidPort_cei[i, j].m_flow) < Modelica.Constants.eps or abs(
        heatPort_cei[i, j].Q_flow) < Modelica.Constants.eps, "(ceiling) there is at least one couple of fluid and thermal connectors coupled");
    end for;
  end for;

  // controls that assure the correctness of the model
  for i in 1:I-1 loop
    for j in 1:J loop
      for k in 1:K loop
        assert(abs(Fe_M[i, j, k] - Fw_M[i + 1, j, k])
           < Modelica.Constants.eps, "Mass flow rates are not equals - x");
      end for;
    end for;
  end for;

  for i in 1:I loop
    for j in 1:J-1 loop
      for k in 1:K loop
        assert(abs(Fn_M[i, j, k] - Fs_M[i, j + 1, k])
           < Modelica.Constants.eps, "Mass flow rates are not equals - y");
      end for;
    end for;
  end for;

  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K - 1 loop
        assert(abs(Fb_M[i, j, k + 1] - Ft_M[i, j, k])
           < Modelica.Constants.eps, "Mass flow rates are not equals - z");
      end for;
    end for;
  end for;

end Room3D_CFDwall_assert;
