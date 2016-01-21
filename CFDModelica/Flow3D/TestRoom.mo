within CFDModelica.Flow3D;
model TestRoom
  extends Modelica.Icons.Example;
  // NUMBER OF VOLUMES
  parameter Integer I=5;
  parameter Integer J=5;
  parameter Integer K=3;
  // SIZES
  parameter Real Width =  10;
  parameter Real Base =   5;
  parameter Real Height = 2.5;
  // WALL TEMPERATURES
  parameter Real T_left = 273.15 + 20;
  parameter Real T_right = 273.15 + 20;
  parameter Real T_floor = 273.15+ 20;
  parameter Real T_cei = 273.15+ 20;
  parameter Real T_front = 273.15+ 20;
  parameter Real T_rear = 273.15+ 20;
  parameter Real Tstart = 273.15 + 20;
  // UNIFORM GRID
  parameter Real x[I]=1/I*ones(I);
  parameter Real y[J]=1/J*ones(J);
  parameter Real z[K]=1/K*ones(K);
  // POWER INPUT
  // (1,3,1) is (left side, middle, bottom)
  parameter Real PowerIN = 100;
  parameter Integer iIN = 1;
  parameter Integer jIN = 3;
  parameter Integer kIN = 1;
  // MEASUREMENT
  // (5,3,2) is (right side, middle, middle)
  parameter Integer iOUT = 5;
  parameter Integer jOUT = 3;
  parameter Integer kOUT = 2;
  CFDModelica.Flow3D.Room3D_CFD room(
    convective=true,
    I=I,
    J=J,
    K=K,
    roomWidth=Width,
    roomBase=Base,
    roomHeight=Height,
    energy=true,
    Tstart=Tstart,
    X_frac=x[:],
    Y_frac=y[:],
    Z_frac=z[:],
    Hright=2,
    Hleft=2,
    Hfloor=2,
    Hcei=2,
    Hfront=2,
    Hrear=2,
    intScheme=CFDModelica.Units.IntScheme.Hybrid,
    intSchemeEnergy=CFDModelica.Units.IntScheme.Hybrid)
    annotation (Placement(transformation(extent={{-42,-62},{18,-2}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow
                                                      PowerSource
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=180,
        origin={38,-110})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature floor(T(displayUnit="K")=
         T_floor)               annotation (Placement(transformation(
          extent={{-100,-130},{-80,-110}})));
  CFDModelica.Thermal.ThermalCollector thermalCollector(
    I=I,
    J=J,
    K=1) annotation (Placement(transformation(extent={{-28,-102},{-8,-82}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature ceiling(T(displayUnit="K")=
         T_cei)                 annotation (Placement(transformation(
          extent={{-100,22},{-80,42}})));
  CFDModelica.Thermal.ThermalCollector cei_coll(
    I=I,
    J=J,
    K=1) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-18,18})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature rear(T(displayUnit="K")=
         T_rear)                annotation (Placement(transformation(
          extent={{-100,58},{-80,78}})));
  CFDModelica.Thermal.ThermalCollector rear_coll(
    I=I,
    J=1,
    K=K) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,68})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature front(T(displayUnit="K")=
         T_front)               annotation (Placement(transformation(
          extent={{-100,-92},{-80,-72}})));
  CFDModelica.Thermal.ThermalCollector front_coll(
    I=I,
    J=1,
    K=K) annotation (Placement(transformation(extent={{-50,-82},{-30,-62}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature left(T(displayUnit="K")=
         T_left)                annotation (Placement(transformation(
          extent={{-100,-36},{-80,-16}})));
  CFDModelica.Thermal.ThermalCollector left_coll1(
    K=K,
    I=1,
    J=J) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-62,-26})));
  CFDModelica.Thermal.ThermalCollector right_coll(
    K=K,
    I=1,
    J=J) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={40,-26})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature right(T(displayUnit="K")=
         T_right)               annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
        rotation=180,
        origin={74,-26})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor temperatureSensor
    annotation (Placement(transformation(extent={{50,-80},{70,-60}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=PowerIN,
    duration=250,
    offset=0,
    startTime=360)
    annotation (Placement(transformation(extent={{90,-120},{70,-100}})));
equation
  // Connections between the input power source, temperature sensor and the room
  connect(PowerSource.port, room.heatPort[iIN, jIN, kIN]) annotation (Line(
      points={{28,-110},{15,-110},{15,-59}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(temperatureSensor.port, room.heatPort[iOUT, jOUT, kOUT]) annotation (Line(
      points={{50,-70},{15,-70},{15,-59}},
      color={191,0,0},
      smooth=Smooth.None));

  // Other connections
  connect(thermalCollector.port_b, floor.port)            annotation (
      Line(
      points={{-18,-102},{-18,-120},{-80,-120}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(thermalCollector.port_a[:, :, 1], room.heatPort_floor)
    annotation (Line(
      points={{-18,-82},{-18,-59}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(ceiling.port, cei_coll.port_b)          annotation (Line(
      points={{-80,32},{-18,32},{-18,28}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(cei_coll.port_a[:, :, 1], room.heatPort_cei)          annotation (
      Line(
      points={{-18,8},{-18,-5}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(rear.port, rear_coll.port_b) annotation (Line(
      points={{-80,68},{-10,68}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(rear_coll.port_a[:, 1, :], room.heatPort_rear) annotation (Line(
      points={{10,68},{10,21},{9,21},{9,-5}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(front_coll.port_b, front.port) annotation (Line(
      points={{-40,-82},{-80,-82}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(front_coll.port_a[:, 1, :], room.heatPort_front) annotation (Line(
      points={{-40,-62},{-40,-53},{-39,-53}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(left.port, left_coll1.port_b) annotation (Line(
      points={{-80,-26},{-76,-26},{-76,-26},{-72,-26}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(left_coll1.port_a[1, :, :], room.heatPort_left) annotation (Line(
      points={{-52,-26},{-39,-26}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(room.heatPort_right, right_coll.port_a[1, :, :]) annotation (Line(
      points={{15,-26},{30,-26}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(right.port, right_coll.port_b) annotation (Line(
      points={{64,-26},{50,-26}},
      color={191,0,0},
      smooth=Smooth.None));

  connect(ramp.y, PowerSource.Q_flow) annotation (Line(
      points={{69,-110},{48,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-150},
            {100,100}}),
                      graphics), Icon(coordinateSystem(
          preserveAspectRatio=true, extent={{-100,-150},{100,100}})),
    experiment(StopTime=3600),
    __Dymola_experimentSetupOutput(equdistant=false, events=false));
end TestRoom;
