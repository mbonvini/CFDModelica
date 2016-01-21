within CFDModelica.Flow3D;
model Room3D_CFD
  extends Room3D_CFDwall_base;
  extends Room3D_CFDwall_init;
  extends Room3D_CFDwall_assert;
  extends Room3D_CFDwall_boundaries;
  extends Room3D_CFDwall_walls;
  extends Room3D_CFDwall_massANDenergy;
  extends Room3D_CFDwall_turbulence;
  extends Room3D_CFDwall_momentum;
  extends Room3D_CFDwall_momentum_boundaries;
equation
  /************************************************************
  **  START OUTPUT EQUATIONS
  **
  ** These are just for post processing purposes
  ************************************************************/
  connect( T_a[:,:,:], T_[:,:,:]);
  connect( rho_a[:,:,:], rho_[:,:,:]);
  connect( Vx_a[:,:,:], Vx_[:,:,:]);
  connect( Vy_a[:,:,:], Vy_[:,:,:]);
  connect( Vz_a[:,:,:], Vz_[:,:,:]);
  connect( P_a[:,:,:],P_[:,:,:]);
  connect( I_a,I_);
  connect( J_a,J_);
  connect( K_a,K_);
  connect( roomBase_a,roomBase_);
  connect( roomWidth_a,roomWidth_);
  connect( roomHeight_a,roomHeight_);
  connect( X_frac_a[:],X_frac_[:]);
  connect( Y_frac_a[:],Y_frac_[:]);
  connect( Z_frac_a[:],Z_frac_[:]);

  for i in 1:I loop
    for j in 1:J loop
      for k in 1:K loop
        T_a[i, j, k] = T[i + 1, j + 1, k +
          1];
        rho_a[i, j, k] = rho[i + 1, j + 1, k
           + 1];
        P_a[i, j, k] = P[i, j, k];
      end for;
    end for;
  end for;
  for i in 1:I+1 loop
    for j in 1:J+2 loop
      for k in 1:K + 2 loop
        Vx_a[i, j, k] = Vx[i, j, k];
      end for;
    end for;
  end for;
  for i in 1:I+2 loop
    for j in 1:J+1 loop
      for k in 1:K + 2 loop
        Vy_a[i, j, k] = Vy[i, j, k];
      end for;
    end for;
  end for;
  for i in 1:I+2 loop
    for j in 1:J+2 loop
      for k in 1:K + 1 loop
        Vz_a[i, j, k] = Vz[i, j, k];
      end for;
    end for;
  end for;
  I_a = I;
  J_a = J;
  K_a = K;
  roomBase_a = roomBase;
  roomHeight_a = roomHeight;
  roomWidth_a  = roomWidth;
  X_frac_a[:] = X_frac[:];
  Y_frac_a[:] = Y_frac[:];
  Z_frac_a[:] = Z_frac[:];
annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-60,40},{40,-60}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-98,12},{78,-2}},
          lineColor={0,0,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Forward,
          textString="%name"),
        Text(
          extent={{-60,-34},{40,-46}},
          lineColor={0,0,0},
          textString="%I x %J x %K"),
        Polygon(
          points={{-60,40},{-40,60},{60,60},{40,40},{-60,40}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{40,-60},{60,-40},{60,60},{40,40},{40,-60}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Line(
          points={{60,80},{40,60}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{84,64},{60,40}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{80,20},{50,20}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{80,-20},{50,-20}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-20,-80},{-20,-60}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{20,-80},{20,-60}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-80,-60},{-46,-26}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-64,-84},{-30,-50}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-80,20},{-60,20}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-80,-20},{-60,-20}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-20,80},{-20,50}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{20,80},{20,50}},
          color={0,0,0},
          smooth=Smooth.None)}));
end Room3D_CFD;
