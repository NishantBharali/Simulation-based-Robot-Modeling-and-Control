close all
clear
clc

% create figure
figure
axis([-4, 4, -4, 4])
grid on
hold on

% save as a video file
v = VideoWriter('Case_2.mp4', 'MPEG-4');
v.FrameRate = 100;
open(v);

% pick your system parameters
m1 = 1;
m2 = 1;
m3 = 1;
I3 = 0.1;
L = 1;
g = 9.81;
deltaT = 0.01;

% initial conditions
theta = [0; 0; 0];
thetadot = [0; 0; 0];
thetadotdot = [0; 0; 0];
time = 0;

% forward kinematics to end-effector
S1 = [0;0;0;1;0;0];
S2 = [0;0;0;0;1;0];
S3 = [0;0;1;0;-L;0];
S = [S1, S2, S3];
M3 = [eye(3), [2*L;0;0]; 0 0 0 1];

% For Case 2
Kp = eye(3)*25;
Kd = eye(3)*25;

M0 = [eye(3), [L;0;0]; 0 0 0 1];
M1 = [eye(3), [L;0;0]; 0 0 0 1];
M2 = [eye(3), [L;0;0]; 0 0 0 1];

for idx = 1:1000

    % get desired position
    theta_d = [-2; 2; pi/4];
    thetadot_d = [0; 0; 0];
    T_d = fk(M3, [S1 S2 S3], theta_d);
    
    % plot the robot
    p0 = [0; 0];
    T1 = fk(M1, S1,theta(1:1,:));
    p1 = T1(1:2,4);                         % position of end of link 1
    T2 = fk(M2,[S1 S2],theta(1:2,:));
    p2 = T2(1:2,4);                         % position of end of link 2
    T3 = fk(M3,[S1 S2 S3],theta(1:3,:));
    p3 = T3(1:2,4);                         % position of end of link 3
    P = [p0, p1, p2, p3];
    cla;
    plot(P(1,:), P(2,:), 'o-', 'color',[1, 0.5, 0],'linewidth',4)
    
    % plot the desired position
    plot(T_d(1,4), T_d(2,4), 'ok', 'MarkerFaceColor','k')
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
    
    % Mass matrix
    M = [m1 + m2 + m3, 0, -L*m3*sin(theta(3));
            0, m2 + m3, L*m3*cos(theta(3));
            -L*m3*sin(theta(3)), L*m3*cos(theta(3)), m3*L^2 + I3];

    % Coriolis matrix
    C = [0, 0, -L*thetadot(3)*m3*cos(theta(3));
        0, 0, -L*thetadot(3)*m3*sin(theta(3));
        0, 0,                 0];

    % Gravity vector
    G = [0; g*m2 + g*m3; L*g*m3*cos(theta(3))];
    
    % Reference from a journal and a book Modern Robotics (Just for reference and knowledge)
    e1 =    theta_d - theta;
    e1dot = thetadot_d - thetadot;

    % Choose your controller tau
    tau = Kp*(theta_d-theta) + Kd*(thetadot_d - thetadot) + G;
    
    % integrate to update velocity and position
    thetadotdot = M \ (tau - C*thetadot - G);
    thetadot = thetadot + deltaT * thetadotdot;
    theta = theta + deltaT * thetadot;
    time = time + deltaT;

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