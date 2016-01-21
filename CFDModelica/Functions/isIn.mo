within CFDModelica.Functions;
function isIn
  "the function returns true if the pair of element (a,b) is contained in M"
  input Integer a "row value";
  input Integer b "column value";
  input Integer[:,2] M = [0,0] "list of positions";
  input Boolean enableALL = false "this flag enable all the connectors";
  output Boolean contained;
protected
  Integer r;
  Integer c;
  Integer i;
algorithm

  contained := false;

  // the sizes of the variable dimension matrix
  r := size(M,1);
  c := size(M,2);

  if enableALL then

    contained := true;

  else

    if r == 0 then
      contained := false;
    else
      i := 1;

      while contained == false and i<=r loop

        if M[i,1] == a then
          // the first element is equal to "a"
          if M[i,2] == b then
            // the second element is equal to "b" -> TRUE
            contained := true;
          else
            // the second element does not match
            i := i+1;
          end if;

        else
          // the first element does not match
          i := i+1;

        end if;

      end while;
    end if;

  end if;
end isIn;
