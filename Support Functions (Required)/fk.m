function F = fk(M_car,S_var,theta_var)

% Forward Kinematics function

    thetas_var = theta_var;
    T1 = eye(4, 4);

    for i=1:length(thetas_var)
       T1 = T1 * expm(bracket_s(S_var(:,i))* thetas_var(i));
    end
    F = T1*M_car;
    function S_matrix = bracket_s(s)
        S_matrix = [0 -s(3) s(2) s(4);
                s(3) 0 -s(1) s(5)
                -s(2) s(1) 0 s(6)
                 0 0 0 0];

    end

end




% function T1 = fk(M_car, S_var, theta_char, L_var)
% 
% omega_var = [0;0;1];
% M_car = [eye(3), [5*L_var;0;0]; 0 0 0 1];
% q1_var = [0;0;0];
% q2_var = [L_var;0;0];
% q3_var = [2*L_var;0;0];
% q4_var = [3*L_var;0;0];
% q5_var = [4*L_var;0;0];
% 
% S1_var = [omega_var; -cross(omega_var, q1_var)];
% S2_var = [omega_var; -cross(omega_var, q2_var)];
% S3_var = [omega_var; -cross(omega_var, q3_var)];
% S4_var = [omega_var; -cross(omega_var, q4_var)];
% S5_var = [omega_var; -cross(omega_var, q5_var)];
% 
% S_var = [S1_var, S2_var, S3_var, S4_var, S5_var];
% 
% T1 = expm(bracket(S1_var) * theta_char(1)) * expm(bracket(S2_var) * theta_char(2)) * expm(bracket(S3_var) * theta_char(3)) * expm(bracket(S4_var) * theta_char(4)) * expm(bracket(S5_var) * theta_char(5)) * M_car;
% 
%     function R_matrix = bracket(S)
%         R_matrix = [ 0 -S(3) S(2) S(4);
%             S(3) 0 -S(1) S(5);
%             -S(2) S(1) 0 S(6);
%             0 0 0 0];
%     end
% end