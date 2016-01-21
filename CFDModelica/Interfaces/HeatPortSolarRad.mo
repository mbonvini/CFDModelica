within CFDModelica.Interfaces;
connector HeatPortSolarRad "Thermal port for 1-dim. heat transfer"
  extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
  annotation(defaultComponentName = "port_a",
    Documentation(info="<HTML>
<p>This connector is used for 1-dimensional heat flow between components.
The variables in the connector are:</p>
<pre>   
   T       Temperature in [Kelvin].
   Q_flow  Heat flow rate in [Watt].
</pre>
<p>According to the Modelica sign convention, a <b>positive</b> heat flow
rate <b>Q_flow</b> is considered to flow <b>into</b> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <b>HeatPort_a</b> and
<b>HeatPort_b</b> are identical with the only exception of the different
<b>icon layout</b>.</p></HTML>
"), Icon(graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid), Rectangle(
          extent={{-60,60},{60,-60}},
          lineColor={0,0,255},
          pattern=LinePattern.None,
          lineThickness=0.5,
          fillColor={255,213,28},
          fillPattern=FillPattern.Solid)}),
    Diagram(graphics));
end HeatPortSolarRad;
