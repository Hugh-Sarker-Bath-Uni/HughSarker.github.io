%% ==============================
%  COUGH DROPLET PLUME MODEL
%  Monte Carlo + turbulence + settling
%% ==============================

clear; clc;

%% ------------------------------
% PARAMETERS
%% ------------------------------
g = 9.81;

density_liquid = 1000;
radius = 5e-6;
air_viscosity = 1.8e-5;

particle_relaxation = (2*density_liquid*radius^2)/(9*air_viscosity);
v_settle = particle_relaxation * g;   % settling speed (m/s)

%% ------------------------------
% MEAN WIND (choose scenario)
%% ------------------------------
U_ref = [3,6,12];      % light, moderate, strong
alpha = [0.13,0.25,0.325]; % field, suburb, urban

scenario_wind = 1;   % moderate
scenario_terrain = 3; % urban

U0 = U_ref(scenario_wind);
a  = alpha(scenario_terrain);

z = 2;
z_ref = 10;

U_mean_scalar = U0 * (z / z_ref)^a;
u_mean = [U_mean_scalar; 0; 0];

%% ------------------------------
% TURBULENCE PARAMETERS
%% ------------------------------
TL = 2;       % turbulence timescale (s)
sigma = 0.4;  % turbulence strength

%% ------------------------------
% TIME SETTINGS
%% ------------------------------
dt = 0.01;
T  = 10;
N  = T/dt;

%% ------------------------------
% MONTE CARLO SETTINGS
%% ------------------------------
Nd = 2000;   % number of droplets

x_end = zeros(3,Nd);

%% ------------------------------
% MAIN SIMULATION LOOP
%% ------------------------------

for i = 1:Nd

    % initial condition
    x = [0;0;2] + 0.02*randn(3,1);

    u_prime = zeros(3,1);

    for n = 1:N-1

        % --------------------------
        % TURBULENCE (OU PROCESS)
        % --------------------------
        xi = randn(3,1);

        u_prime = u_prime - (dt/TL)*u_prime + sigma*sqrt(dt)*xi;

        % --------------------------
        % MEAN FLOW + TURBULENCE
        % --------------------------
        u_air = u_mean + u_prime;

        dxdt = u_air;

        % settling (downwards)
        dxdt(3) = dxdt(3) - v_settle;

        % --------------------------
        % UPDATE POSITION
        % --------------------------
        x_new = x + dxdt*dt;

        % --------------------------
        % DEPOSITION CONDITION
        % --------------------------
        if x_new(3) <= 0
            x_new(3) = 0;   % lock to ground
            x = x_new;
            break;          % stop this droplet
        end

        x = x_new;

    end

    x_end(:,i) = x;

end

%% ==============================
% PLOTTING: PLUME DISTRIBUTION
%% ==============================

%% 1. Ground footprint (XY plane)
figure;
scatter(x_end(1,:), x_end(2,:), 8, 'filled');

xlabel('x (m)');
ylabel('y (m)');
title('Cough Plume Footprint (Ground Level)');
grid on;
axis equal;

%% 2. 3D plume
figure;
scatter3( ...
    x_end(1,:), ...
    x_end(2,:), ...
    x_end(3,:), ...
    10, ...
    x_end(3,:), ...
    'filled' ...
);

colorbar;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('3D Cough Droplet Plume (Monte Carlo)');
grid on;
view(3);

%% 3. Density map (XY histogram)
figure;
hist3([x_end(1,:)' x_end(2,:)'], 'Nbins',[60 60]);

xlabel('x (m)');
ylabel('y (m)');
zlabel('Density');
title('Plume Density Field (XY Projection)');


%what we see here is that the mean wind speed pushed the particle droplets along the x axis (directionally away from the mouth). Turbulence causes the droplet to stray away the further away from the inital position given that it is moving upwards as if it heads downwards it will stop at the ground. this is why you see a sharp rise upwards as we move along the x axis.

x_mean = mean(x_end, 2);

x_mean


%x: from this we can infer that the wind is strang and very affective at carrying the droplet away from the mouth


%y: this answer is approx zero, which makes sense as turbulence has symettric randomness due the the norm dist.

%z: this is a suprising result because it shows that on average the droplets end up higher

ground_hits = x_end(3,:) == 0;

x_ground_mean = mean(x_end(1,ground_hits));
y_ground_mean = mean(x_end(2,ground_hits));

scatter(x_end(1,:), x_end(2,:), 10, 'filled');
title('Droplet landing')
mean_ground = [x_ground_mean, y_ground_mean]


Nd = 50;              % number of droplets
x_all = zeros(3,N,Nd);

for i = 1:Nd

    x = zeros(3,N);
    x(:,1) = [0;0;2]; % start at a height of 2m (very tall person)

    u_prime = zeros(3,1);

    for n = 1:N-1

        xi = randn(3,1);

        u_prime = u_prime ...
            - (dt/TL)*u_prime ...
            + sigma*sqrt(dt)*xi;

        z_eff = max(x(3,n),0.01);

        U_mean = U0*(z_eff/z_ref)^a;
        u_mean = [U_mean;0;0];

        u_air = u_mean + u_prime;

        dxdt = u_air;
        dxdt(3) = dxdt(3) - v_settle;

        x(:,n+1) = x(:,n) + dxdt*dt;

        if x(3,n+1) <= 0
            x(3,n+1:end) = 0;
            break
        end

    end

    x_all(:,:,i) = x;

end

figure;
hold on
grid on

for i = 1:Nd
    plot3(x_all(1,:,i), x_all(2,:,i), x_all(3,:,i))
end

xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Multiple Cough Droplet Trajectories (Urban-Moderate-Wind Model)');
view(1);

hold on





