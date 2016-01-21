within CFDModelica.Units;
type Method = enumeration(
    Diff,
    FV,
    FVpatankar,
    FVbook) "Method employed for the discretisation";
