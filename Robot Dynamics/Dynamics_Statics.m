clc
clear

syms theta1 theta2 theta3 real

% What wrench does the robot need to apply at the end-effector to maintain static equilibrium?

% Solution: (Check the PDF demonstration for better visualization of the outputs)

% Opposite force of +15 N over positive y-axis
% Thus the wrench required is provided below with the formulae: Fb = [mb;fb] and moment is zero as the force is directly applied over the body frame

Fb = [0;0;0;0;15;0];

% Case 1: Let 𝐿 = 1 and let 𝜃 = [0, 𝜋/4, 𝜋/4] ' . Find the joint torques needed to balance out the force applied by the human. 
% Case 2: Let 𝐿 = 1 and let 𝜃 = [0, 𝜋/8, 0] ' . Find the joint torques needed to balance out the force applied by the human

theta = [theta1;theta2;theta3];


omega = [0;0;1];

q1 = [0;0;0];
q2 = [1;0;0];
q3 = [2;0;0];

S1 = [omega; -cross(omega, q1)];
S2 = [omega; -cross(omega, q2)];
S3 = [omega; -cross(omega, q3)];

S_eq = [S1, S2, S3];   
M = [eye(3), [3;0;0]; 0 0 0 1];

% T with initial joint positions
T_0 = simplify(expand(fk(M, S_eq, theta)))
R_0 = T_0(1:3, 1:3);
JS = simplify(expand(JacS(S_eq, theta))) %Space Jacobian
Jb = simplify(expand(adjointM(inv(T_0))*JS)) %Body Jacobian
J_geometric = simplify(expand([R_0, zeros(3); zeros(3), R_0] * Jb)) %Geometric Jacobian

% Fs calculation: 
Fs = simplify(expand(adjointM(inv(T_0)))'*Fb)

% symbolic tau calculation
tau = simplify(expand(JS'*Fs))

% CASE 1
Case1_tau = double(subs(tau,[theta1,theta2,theta3],[0,pi/4,pi/4]))

% CASE 2
Case2_tau = double(subs(tau,[theta1,theta2,theta3],[0,pi/8,0]))

% Let ∥𝜏∥ be the magnitude (i.e., the length) of the joint torque vector. 
%     • Find a joint position 𝜃 that maximizes ∥𝜏∥ 
%     • Find a joint position 𝜃 that minimizes ∥𝜏∥ 

% ||tau||
magnitude_tau = simplify(expand(norm(tau)))


% Maximum ||tau|| with theta1 = theta2 = theta3 = 0
% or theta1 = 90 degrees, theta2 = 0, theta3 = 0
Max__mag_tau = double(subs(magnitude_tau,[theta1,theta2,theta3],[0,0,0]))
% f = (15*sqrt(sym(2))*sqrt(cos(2*theta3) + 4*cos(theta3) + 2*(cos(theta2 + theta3) + cos(theta3) + 1)^2 + 5))/2

% Minimum ||tau|| with theta1 = 90 degrees, theta2 = 90 degrees, theta3 = 180 degrees
Min_mag_tau = double(subs(magnitude_tau,[theta1,theta2,theta3],[-pi/2,pi/2,pi]))

% Verifying the optimal results of all theta values for maximum and minimum
% evaluation of the magnitude of joint torque values:

fprintf(['Verifying results for maximum and minimum magnitude of tau with' ...
    'optimal theta values:\n']);

% Define the objective function to maximize Magnitude_tau
objectiveFunction_Max = @(theta) -double(norm(subs(magnitude_tau, [theta1, theta2, theta3], double(theta))));
% Define the objective function to mainimize Magnitude_tau
objectiveFunction_Min = @(theta) double(norm(subs(magnitude_tau, [theta1, theta2, theta3], double(theta))));

% Define initial guess for thetas
x0 = [pi/4, pi/4, pi/4];

% Define bounds on thetas
lb = [0, 0, 0];
ub = [pi, pi, pi];

% Set up the optimization options
options = optimoptions('fmincon', 'Display', 'iter'); % Display optimization process

% Solve the optimization problem to find maximum magnitude_tau
[xMax, fMax] = fmincon(objectiveFunction_Max, x0, [], [], [], [], lb, ub, [], options);

% Display the results (Maximum)
fprintf('Maximum Magnitude_tau: %f\n', abs(-fMax)); % Negate the value back to the original form
fprintf('Optimal values for theta1, theta2, and theta3: %f, %f, %f\n', rad2deg(xMax(1)), rad2deg(xMax(2)), rad2deg(xMax(3)));

% Solve the optimization problem to find minimum magnitude_tau
[xMin, fMin] = fmincon(objectiveFunction_Min, x0, [], [], [], [], lb, ub, [], options);

% Display the results (Minimum)
fprintf('Minimum Magnitude_tau: %f\n', abs(-fMin)); % Negate the value back to the original form
fprintf('Optimal values for theta1, theta2, and theta3: %f, %f, %f\n', rad2deg(xMin(1)), rad2deg(xMin(2)), rad2deg(xMin(3)));


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