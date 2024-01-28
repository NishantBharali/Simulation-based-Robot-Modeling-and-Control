clc;
clear;

syms L g m1 m2 m3 theta1 theta2 theta3 theta1_dot theta2_dot theta3_dot theta1_dot_dot theta2_dot_dot theta3_dot_dot Ix1 Ix2 Ix3 Iy1 Iy2 Iy3 Iz1 Iz2 Iz3 real

theta = [theta1; theta2; theta3];
thetadot = [theta1_dot; theta2_dot; theta3_dot];
thetadotdot = [theta1_dot_dot; theta2_dot_dot; theta3_dot_dot];

% Home matrix till m1 center of mass
M1 = [eye(3),[L;0;0]; 0 0 0 1];
% Home matrix till m2 center of mass
M2 = [eye(3), [L;0;0]; 0 0 0 1];
% Home matrix till m3 center of mass
M3 = [eye(3), [2*L;0;0]; 0 0 0 1];

% Screw for joint 1
S1 = [0;0;0;1;0;0];
% Screw for joint 2
S2 = [0;0;0;0;1;0];
% Screw for joint 3
S3 = [0;0;1;0;0;0];

% Considering m1 center of mass as end-effector
S_eq1 = [S1, zeros(6, 1), zeros(6, 1)];
% Considering m2 center of mass as end-effector
S_eq2 = [S1, S2, zeros(6, 1)];
% Considering m3 center of mass as end-effector
S_eq3 = [S1, S2, S3];

% For center of mass m1
% Forward kinematics
T_1 = fk(M1, S_eq1, theta)
R_1 = T_1(1:3, 1:3);
% Space Jacobian
Js_1 = simplify(expand(JacS(S_eq1, theta)));
% Body Jacobian
Jb_1 = adjointM(inv(T_1))*Js_1;
% Geometric Jacobian
J_geometric_1 = simplify(expand([R_1, zeros(3); zeros(3), R_1] * Jb_1));
% NOTE: For Jw(x1:y1, x2:y2) and Jv(x1:y1, x2:y2), number of columns y1 and
% y2 vary according to the number of joints, thus we change accordingly
Jw1 = J_geometric_1(1:3,1:3)
Jv1 = J_geometric_1(4:6, 1:3)
Inertia_1 = [[Ix1 0 0]
    [0 Iy1 0]
    [0 0 Iz1]];

% For m2 center of mass
T_2 = fk(M2, S_eq2, theta)
R_2 = T_2(1:3, 1:3);
Js_2 = simplify(expand(JacS(S_eq2, theta))); % Space Jacobian
Jb_2 = adjointM(inv(T_2))*Js_2; % Body Jacobian
J_geometric_2 = simplify(expand([R_2, zeros(3); zeros(3), R_2] * Jb_2)); % Geometric Jacobian
Jw2 = J_geometric_2(1:3,1:3)
Jv2 = J_geometric_2(4:6, 1:3)

Inertia_2 = [[Ix2 0 0]
    [0 Iy2 0]
    [0 0 Iz2]];

% For m3 center of mass
T_3 = fk(M3, S_eq3, theta)
R_3 = T_3(1:3, 1:3);
Js_3 = simplify(expand(JacS(S_eq3, theta))); % Space Jacobian
Jb_3 = adjointM(inv(T_3))*Js_3; % Body Jacobian
J_geometric_3 = simplify(expand([R_3, zeros(3); zeros(3), R_3] * Jb_3)); % Geometric Jacobian
Jw3 = J_geometric_3(1:3,1:3)
Jv3 = J_geometric_3(4:6, 1:3)

Inertia_3 = [[Ix3 0 0]
    [0 Iy3 0]
    [0 0 Iz3]];


% Mass Matrix
Mass_Matrix = simplify(expand(m1*(Jv1'*Jv1) + Jw1'*R_1*Inertia_1*R_1'*Jw1 + m2*(Jv2'*Jv2) + Jw2'*R_2*Inertia_2*R_2'*Jw2 + m3*(Jv3'*Jv3) + Jw3'*R_3*Inertia_3*R_3'*Jw3))

% Coriolis Matrix
Coriolis_Matrix = coriolis(Mass_Matrix, theta, thetadot)

% Height of each center of mass
h1 = T_1(2, 4)
h2 = T_2(2, 4)
h3 = T_3(2, 4)

% Potential Energy
P = g*m1*h1 + g*m2*h2 + g*m3*h3

% Gravity vector
gravity_vector = simplify(expand([diff(P, theta(1)); diff(P, theta(2)); diff(P, theta(3))]))

% Tau calculation
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