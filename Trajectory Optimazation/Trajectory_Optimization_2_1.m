clear
clc
close all

% Start and Goal orientations
theta_start = [0;0];
theta_goal = [1;1];

% Initial trajectory variables

n = 2; % No. of joints / 2-D trajectory. 
k = 10; % No. of waypoints

% Obstacle 1 Paramters
r1 = 0.3;
center1 = [0.55;0.5];

xi_0 = zeros(n, k);
xi_0_vec = reshape(xi_0, [],1);

% Equality constraints for start and goal positions
A = [eye(n), zeros(n,n*(k-1)) ;...
    zeros(n,n*(k-1)), eye(n) ];

B = [theta_start;theta_goal];

% Nonlinear optimization
options = optimoptions('fmincon','Display','iter',...   
    'Algorithm','sqp','MaxFunctionEvaluations',1e4);
xi_star_vec = fmincon(@(xi) cost(xi), xi_0_vec, ...
    [], [], A, B, [], [], [], options);

xi_star = reshape(xi_star_vec,2,[]); % to implement

figure
grid on
hold on
axis([0, 1, 0, 1])
axis equal
viscircles(center1', r1, 'Color', [0.5, 0.5, 0.5]);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k')
plot(1, 1, 'ko', 'MarkerFaceColor', 'k')

% Plot Result
grid on
hold on
axis equal
plot(xi_star(1,:), xi_star(2,:), 'o-',...
    'Color', [1, 0.5, 0], 'LineWidth', 3);

% Cost function to minimize

function C = cost(xi)
        
    gamma = 20;
    xi = reshape(xi,2,[]);
    C = 0;
    r1 = 0.3;
    center1 = [0.55;0.5];
    
    for idx = 2:length(xi)
        Urep1 = 0;

        if (norm(center1 - xi(:,idx)) <=r1) 
            Urep1 = 0.5*gamma*((1/(norm(center1-xi(:,idx)))) - (1/r1))^2;
        end

        C = Urep1  + C + norm(xi(:,idx) - xi(:,idx-1))^2;
    end
end