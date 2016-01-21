within CFDModelica.Interfaces;
connector flangeB "Fluid port for fluid transfer"
  extends Modelica.Fluid.Interfaces.FluidPort;
  annotation (Icon(graphics={Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={85,170,255},
          fillPattern=FillPattern.Solid), Ellipse(
          extent={{-62,60},{60,-58}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid)}));
end flangeB;
