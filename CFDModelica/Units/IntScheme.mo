within CFDModelica.Units;
type IntScheme = enumeration(
    Upwind,
    CD,
    Hybrid,
    PowerLaw,
    QUICK) "Method employed for the integration of terms";
