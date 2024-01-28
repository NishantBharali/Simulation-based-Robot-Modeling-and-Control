    function Adj1 = adjointM(T1)

    %Adjoint function

        R = T1(1:3, 1:3);
        p = T1(1:3, 4);
        Adj1 = [R, zeros(3,3); bracket3(p)*R, R];
        
    end