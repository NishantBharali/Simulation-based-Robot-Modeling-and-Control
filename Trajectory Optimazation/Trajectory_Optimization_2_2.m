
clear
clc
close all

% Start and Goal orientations
theta_start = [0;0];
theta_goal = [1;1];

% Initial trajectory variables
n = 2; % No. of joints/ 2-D trajectory
k = 15; % No. of waypoints

% Obstacles 1 and 2 parameters
% First obstacle's radius and center
r1 = 0.3;
center1 = [0.5;0.3];

% Second obstacle's radius and center
r2 = 0.2; 
center2 = [0.5;0.7];

xi_0 = zeros(n, k); % Initial trajectory

xi_0_vec = reshape(xi_0, [],1); % Reshape for the need of optimization

% Equality constraints for start and goal positions
A = [eye(n), zeros(n,n*(k-1)) ;...
 zeros(n,n*(k-1)), eye(n) ];
B = [theta_start;theta_goal];

% Nonlinear optimization
options = optimoptions('fmincon','Display','iter',...
 'Algorithm','sqp','MaxFunctionEvaluations',1e5);

xi_star_vec = fmincon(@(xi) cost(xi), xi_0_vec, ...
 [], [], A, B, [], [], [], options);
xi_star = reshape(xi_star_vec,2,[]); % final optimized trajectory

% Plot obstacles
figure
grid on
hold on
axis([0, 1, 0, 1])
axis equal
viscircles(center1', r1, 'Color', [0.5, 0.5, 0.5]);
viscircles(center2', r2, 'Color', [0.5, 0.5, 0.5]);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k')
plot(1, 1, 'ko', 'MarkerFaceColor', 'k')
% Plot result
grid on
hold on
axis equal
plot(xi_0(1,:), xi_0(2,:), 'o-', 'Color', [0.3, 0.3, ...
 0.3], 'LineWidth', 3);
plot(xi_star(1,:), xi_star(2,:), 'o-',...
 'Color', [1, 0.5, 0], 'LineWidth', 3);

% Cost function to minimize

function C = cost(xi)
 gamma=20;
 xi = reshape(xi,2,[]);
 C = 0;
 r1 = 0.3;
 center1 = [0.5;0.3];
 r2 = 0.2;
 center2 = [0.5;0.7];
 Urep1 = 0;
 Urep2 = 0;
 for idx = 2:length(xi)
 
 % First obstacle
 if (norm(center1 - xi(:,idx)) <= r1) 
 Urep1 = 0.5*gamma*((1/(norm(center1 - xi(:,idx)))) - (1/r1))^2;
 end
 % Second obstacle
 if (norm(center2 - xi(:,idx)) <= r2) 
 Urep2 = 0.5*gamma*((1/(norm(center2 - xi(:,idx)))) - (1/r2))^2;
 end

 % Total cost
 C = C + norm(xi(:,idx) - xi(:,idx-1))^2 + Urep1 + Urep2;
 end
end