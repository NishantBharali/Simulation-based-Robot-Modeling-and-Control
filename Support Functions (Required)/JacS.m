function  Js = JacS(S,theta)
    
% Main Space Jacobian Function

    T = eye(4);
    Js=sym(zeros(6, length(theta)));

    for i = 1:length(theta)
        Si = S(:,i);
        Js(:,i) = Adjoint(T) * Si;
        T = T * expm(bracket_s(Si) * theta(i));

    end

    function S_matrix = bracket_s(s)
        S_matrix = [0 -s(3) s(2) s(4);
                s(3) 0 -s(1) s(5)
                -s(2) s(1) 0 s(6)
                 0 0 0 0];

    end

    function Adt = Adjoint(T)
        
        R = T(1:3,1:3);
        p = T(1:3,4);
        Adt = [R,zeros(3,3);bracket(p)*R,R];

        function w_matrix= bracket(w)
            w_matrix = [0 -w(3) w(2);
                        w(3) 0 -w(1);
                        -w(2) w(1) 0];

        end
    end
    
end


