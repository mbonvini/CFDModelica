within CFDModelica.Functions;
function A
  input Real Peclet;
  input CFDModelica.Units.IntScheme s;
  output Real A;
algorithm
  if s == CFDModelica.Units.IntScheme.Upwind then
    A := 1;
  elseif s == CFDModelica.Units.IntScheme.CD then
    A := 1 - 0.5*abs(Peclet);
  elseif s == CFDModelica.Units.IntScheme.Hybrid then
    A := max(0,1-0.5*abs(Peclet));
  elseif s == CFDModelica.Units.IntScheme.PowerLaw then
    A := max(0,(1-0.5*abs(Peclet))^5);
  elseif s == CFDModelica.Units.IntScheme.QUICK then
    A := 1;
  else
    A:= 1;
  end if;
end A;
