close all
clear
clc

% create figure
figure
axis([-6, 6, -6, 6])
grid on
hold on

% save as a video file
v = VideoWriter('Case_1.mp4', 'MPEG-4');
v.FrameRate = 25;
open(v);

%initial joint values
L = 1;
theta = [pi/8; pi/8; pi/8; pi/8; pi/8];

omega = [0;0;1];

S1 = [0 0 1 0 0 0]';
S2 = [0 0 1 0 -1*L 0]';
S3 = [0 0 1 0 -2*L 0]';
S4 = [0 0 1 0 -3*L 0]';
S5 = [0 0 1 0 -4*L 0]';

S_eq = [S1, S2, S3, S4, S5];   
M = [eye(3), [5*L;0;0]; 0 0 0 1];
M1 = [eye(3), [1*L;0;0]; 0 0 0 1];
M2 = [eye(3), [2*L;0;0]; 0 0 0 1];
M3 = [eye(3), [3*L;0;0]; 0 0 0 1];
M4 = [eye(3), [4*L;0;0]; 0 0 0 1];

% Given desired Transformation matrices T_d
T_d = [rotz(pi/4), [3;2;0]; 0 0 0 1];
Xd = [r2axisangle(T_d(1:3, 1:3)); T_d(1:3,4)];

% T with initial joint positions
T = fk(M, S_eq, theta);
X = [r2axisangle(T(1:3, 1:3)); T(1:3,4)];

while norm(Xd - X) > 1e-2
% plot the robot
% 1. get the position of each link
    p0 = [0; 0];
    T1 = fk(M1, S1, theta(1));
    T2 = fk(M2, [S1, S2], [theta(1), theta(2)]);
    T3 = fk(M3, [S1, S2, S3], [theta(1), theta(2), theta(3)]);
    T4 = fk(M4, [S1, S2, S3, S4], [theta(1), theta(2), theta(3), theta(4)]);
    P_v = [p0, T1(1:2, 4), T2(1:2, 4), T3(1:2, 4), T4(1:2, 4), T(1:2, 4)];

% 2. draw the robot and save the frame

    cla;
    plot(P_v(1,:), P_v(2,:), 'o-', 'color',[1, 0.5, 0],'linewidth',4)
    drawnow
    frame = getframe(gcf);
    writeVideo(v, frame);

% My code Implementation
    JS = JacS(S_eq, theta); % Updating Space Jacobian
    Jb = adjointM(inv(T))*JS; % Updating Body Jacobian
    J_geometric = [T(1:3, 1:3) zeros(3); zeros(3) T(1:3, 1:3)] * Jb; % Updated Geometric Jacobian
    V = Xd - X;

    % Here, we set b vector as the following: b = [-theta(1);0;0;0;0]
    delta_theta = pinv(J_geometric)*V;

    % +(eye(5) - pinv(J_geometric)*J_geometric)*[0;0;0;0;0];
    
    % Updating theta until the while loop is satisfied to get the desired inverse kinematics (joint positions), thus simulating the robot
    theta = double(theta + 0.1 * delta_theta);
    T = fk(M, S_eq, theta);
    X = [r2axisangle(T(1:3, 1:3)); T(1:3,4)];

end

close(v);
close all

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