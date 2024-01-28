clc;
clear;

syms L g m1 m2 theta1 theta2 theta1_dot theta2_dot theta1_dot_dot theta2_dot_dot Ix1 Ix2 Ix3 Iy1 Iy2 Iy3 Iz1 Iz2 real

theta = [theta1; theta2];
thetadot = [theta1_dot; theta2_dot];
thetadotdot = [theta1_dot_dot; theta2_dot_dot];

% home matrix for center of mass m1
M1 = [roty(0),[0;L;-L/2]; 0 0 0 1]; 

% home matrix for center of mass m2
M2 = [eye(3), [0;L;-L]; 0 0 0 1];

S1 = [0;1;0;0;0;0];
S2 = [0;0;0;0;0;-1];

S_eq1 = [S1, [0;0;0;0;0;0]];
S_eq2 = [S1, S2];

% For center of mass m1
T_1 = fk(M1, S_eq1, theta)
R_1 = T_1(1:3, 1:3);
JS_1 = simplify(expand(JacS(S_eq1, theta))); % Space Jacobian
Jb_1 = adjointM(inv(T_1))*JS_1; % Body Jacobian
J_geometric_1 = simplify(expand([R_1, zeros(3); zeros(3), R_1] * Jb_1)); % Geometric Jacobian
Jw1 = J_geometric_1(1:3,1:2)
Jv1 = J_geometric_1(4:6, 1:2)
Inertia_1 = [[Ix1 0 0]
    [0 Iy1 0]
    [0 0 Iz1]];


T_2 = fk(M2, S_eq2, theta)
R_2 = T_2(1:3, 1:3);
JS_2 = simplify(expand(JacS(S_eq2, theta))); % Space Jacobian
Jb_2 = adjointM(inv(T_2))*JS_2; % Body Jacobian
J_geometric_2 = simplify(expand([R_2, zeros(3); zeros(3), R_2] * Jb_2)); % Geometric Jacobian
Jw2 = J_geometric_2(1:3,1:2)
Jv2 = J_geometric_2(4:6, 1:2)

Inertia_2 = [[Ix2 0 0]
    [0 Iy2 0]
    [0 0 Iz2]];

% Mass matrix evaluation
Mass_Matrix = simplify(expand(m1*(Jv1'*Jv1) + Jw1'*R_1*Inertia_1*R_1'*Jw1 + m2*(Jv2'*Jv2) + Jw2'*R_2*Inertia_2*R_2'*Jw2))

% Coriolis matrix evaluation
Coriolis_Matrix = coriolis(Mass_Matrix, theta, thetadot)

% Height evaluations as it is acting along the z-axis
h1 = T_1(3, 4);
h2 = T_2(3, 4);

% Potential energy evaluation
P = g*m1*h1 + g*m2*h2;

% Gravity vector evaluation
gravity_vector = simplify(expand([diff(P, theta1); diff(P, theta2)]))

% Final tau (joint torques) evaluation
tau = simplify(expand(Mass_Matrix*thetadotdot + Coriolis_Matrix*thetadot + gravity_vector))


% Supporting Functions
function R_matrix = bracket3(S)
    R_matrix = [0 -S(3) S(2); S(3) 0 -S(1); -S(2) S(1) 0];
end

function R_matrix = bracket4(S)
            R_matrix = [ 0 -S(3) S(2) S(4);
            S(3) 0 -S(1) S(5);
            -S(2) S(1) 0 S(6);
            0 0 0 0];
end

function F = fk(M_car,S_var,theta_var)

% Forward Kinematics function

    thetas_var = theta_var;
    T1 = eye(4, 4);

    for i=1:length(thetas_var)
       T1 = T1 * expm(bracket4(S_var(:,i))* thetas_var(i));
    end
    F = T1*M_car;

end

function Adj1 = adjointM(T1)

%Adjoint function

    R = T1(1:3, 1:3);
    p = T1(1:3, 4);
    Adj1 = [R, zeros(3,3); bracket3(p)*R, R];
        
end


function omega_r2 = r2axisangle(R)
    if norm(R - eye(3)) < 1e-3
        omega_r2 = [0;0;0];
    else
        theta_axis = acos(0.5 * (trace(R) - 1));
        omega_hat = 1/(2*sin(theta_axis)) * (R - R');
        omega_hat = [omega_hat(3,2);
        omega_hat(1,3);
        omega_hat(2,1)];
        omega_r2 = omega_hat * theta_axis;
    end
end

function  Js = JacS(S,theta)
    
% Main Space Jacobian Function

    T = eye(4);
    Js=sym(zeros(6, length(theta)));

    for i = 1:length(theta)
        Si = S(:,i);
        Js(:,i) = adjointM(T) * Si;
        T = T * expm(bracket4(Si) * theta(i));

    end
    
end

function cmatrix = coriolis(m, theta, thetadot)

    n = length(theta); % Depends upon the no. of joints
    cmatrix = sym(zeros(size(m))); % Pre-allocating and initializing the matrix


    for k = (1:(size(cmatrix,1)))
        sum = 0; % Initializing the coriolis
        for j = (1:(size(cmatrix,2)))
            for i = (1:n)
                sum = sum + 1/2*(gradient(m(k,j),theta(i)) + gradient(m(k,i),theta(j)) - gradient(m(i,j),theta(k)))*thetadot(i);
            end
            cmatrix(k, j) = sum;
            sum = 0;
        end
    end
end
