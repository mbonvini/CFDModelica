within CFDModelica.Interfaces;
connector flangeA "Fluid port for fluid transfer"
  extends Modelica.Fluid.Interfaces.FluidPort;
  annotation (Icon(graphics={Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={85,170,255},
          fillPattern=FillPattern.Solid)}));
end flangeA;
