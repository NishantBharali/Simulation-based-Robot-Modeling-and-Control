clear
close all

% Start and goal environments
theta_start= [0; 0];
theta_goal = [1; 1];

% Obstacle parameters
obs_c21= [0.2; 0.5];
obs_r21= 0.125;
obs_c22= [0.7; 0.5];
obs_r22= 0.225;

% Visualizing the environment
figure
grid on
hold on
axis([0, 1, 0, 1])
axis equal
viscircles(obs_c21', obs_r21, 'Color', [0.5, 0.5, 0.5]);
viscircles(obs_c22', obs_r22, 'Color', [0.5, 0.5, 0.5]);
plot(0, 0, 'ko', 'MarkerFaceColor', 'k')
plot(1, 1, 'ko', 'MarkerFaceColor', 'k')

%variables
alpha= 0.01;
epsilon = 0.1;
delta= 0.01;

% initial trajectory
theta(:,1) = theta_start;
t=1;
del_Unet=1;

while norm(del_Unet)> epsilon
    del_Ux= U_theta(theta(:,t) + [delta;0]);
    del_Uy= U_theta(theta(:,t)+ [0;delta]);
    del_U= U_theta(theta(:,t));

    del_Unet=[del_Ux-del_U; del_Uy-del_U]/delta;

    theta(:,t+1)= theta(:,t)- alpha*del_Unet;
    t=t+1;
end


grid on
hold on
axis equal
plot(theta(1,:), theta(2,:), 'o-',...
    'Color', [1, 0.5, 0], 'LineWidth', 2);



function U = U_theta(theta) 
    beta=2;
    gamma=1;
    theta_goal = [1; 1];
    obs_c21= [0.2; 0.5];
    obs_r21= 0.125;
    obs_c22= [0.7; 0.5];
    obs_r22= 0.225;

    Urep1=0;
    Urep2=0;


    Uatt=0.5*beta*norm(theta_goal-theta)^2;    

    if norm(obs_c21-theta)<= obs_r21
        Urep1=0.5* gamma*((1/norm(obs_c21-theta))- (1/obs_r21))^2;
    end
    if norm(obs_c22-theta)<= obs_r22
        Urep2= 0.5*gamma*((1/norm(obs_c22-theta))- (1/obs_r22))^2;
    end

    Urep= Urep1+ Urep2;
    disp(Urep)
    U= Uatt + Urep;
end