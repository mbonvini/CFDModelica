within CFDModelica.Thermal;
model ThermalCollector "Collects IxJxK heat flows"
  parameter Integer I(min=1)=3 "Number of collected heat flows";
  parameter Integer J(min=1)=3 "Number of collected heat flows";
  parameter Integer K(min=1)=3 "Number of collected heat flows";
  CFDModelica.Interfaces.HeatPort port_a[I,J,K]
    annotation (Placement(transformation(extent={{-10,110},{10,90}})));
  CFDModelica.Interfaces.HeatPort port_b
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));

equation
  port_b.Q_flow + sum(port_a.Q_flow) = 0;
  port_a.T = fill(port_b.T, I, J, K);
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}}), graphics),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Text(
          extent={{-150,-30},{150,-70}},
          textString="%name",
          lineColor={0,0,255}),
        Text(
          extent={{-142,70},{148,52}},
          lineColor={0,0,0},
          textString="I=%I, J=%J, K=%K"),
        Line(
          points={{0,90},{0,40}},
          color={181,0,0},
          smooth=Smooth.None),
        Rectangle(
          extent={{-60,40},{60,30}},
          lineColor={181,0,0},
          fillColor={181,0,0},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-60,30},{0,-30},{0,-90}},
          color={181,0,0},
          smooth=Smooth.None),
        Line(
          points={{0,-30},{-20,30}},
          color={181,0,0},
          smooth=Smooth.None),
        Line(
          points={{0,-30},{20,30}},
          color={181,0,0},
          smooth=Smooth.None),
        Line(
          points={{0,-30},{60,30}},
          color={181,0,0},
          smooth=Smooth.None)}),
    Documentation(info="<html>
<p>
This is a model to collect the heat flows from <i>m</i> heatports to one single heatport.
</p>
</html>"));
end ThermalCollector;
