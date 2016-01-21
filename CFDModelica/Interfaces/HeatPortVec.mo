within CFDModelica.Interfaces;
connector HeatPortVec "Thermal port for 1-dim. heat transfer"
  parameter Integer n "lumps";
  Modelica.SIunits.Temperature T[n] "Port temperature";
  flow Modelica.SIunits.HeatFlowRate Q_flow[n]
    "Heat flow rate (positive if flowing from outside into the component)";
  annotation (Icon(graphics={Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid), Text(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          textString=
               "%n")}));
end HeatPortVec;
