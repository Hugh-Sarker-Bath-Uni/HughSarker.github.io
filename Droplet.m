%modelling droplets

% we begin with the small droplet. first we consider the particle ralexation timescale:
%-------------------------------------------------------------------------
g = 9.81
density_liquid = 1000
radius = 5 * 10^-6
air_viscosity = 1.8 * 10^-5

particle_relaxation = (2*density_liquid*radius^2)/(9*air_viscosity)

%-------------------------------------------------------------------------

%partical relaxation is proven to be tiny so we assume velocity of the particle is equal to  u_air

%now to account for the setting velocity

v_settle = -particle_relaxation*g


% i know choose to model the air flow based on mean background airflow


%mean Air flow
U_ref = [3,6,12] %light , moderate, heavy
z = 2
z_ref = 10
alpha = [0.13, 0.25, 0.325] %field, suburb, urban

U0 = U_ref(1)
a = alpha(3)


%-------------------------------------------------------------------------

% turbulence
TL = 2;          % turbulence timescale
sigma = 0.4;     % intensity
dt = 0.01;
T = 10;
N = T/dt;

u_prime = zeros(3,N);
u_air = zeros(3,N);
x = zeros(3,N);
x(:,1) = [0;0;2]; % initial position (2 m height)
U_mean = U0 * (x(3,1) / z_ref)^a;
u_mean = [U_mean; 0; 0];


for n = 1:N-1

    % --- turbulence update ---
    xi = randn(3,1); %randomness added but based on the standard normal distn this is how turbulence should be modelled

    u_prime(:,n+1) = u_prime(:,n) ...
        - (dt/TL)*u_prime(:,n) ...
        + sigma*sqrt(dt)*xi;

    % --- total wind ---
    u_air(:,n) = u_mean + u_prime(:,n);

end

%-------------------------------------------------------------------------

for n = 1:N-1

    % velocity = airflow + settling
    dxdt = u_air(:,n);

    % settling only in z-direction (downwards)
    dxdt(3) = dxdt(3) - v_settle;

    % Euler update
    x(:,n+1) = x(:,n) + dxdt*dt;

end

%so what we have done here is that we have create a ODE that says the velocity is based on the mean air flow plus some turbulence, which is random leading us to have different graphs each time.




%----------------------------------------
    %the plot of trajectoy%
%----------------------------------------

figure;
plot3(x(1,:), x(2,:), x(3,:), 'b-', 'LineWidth', 1.5);
grid on;

xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Cough Droplet Trajectory (Urban Model)');
view(1);

%----------------------------------------
    %the plot of height%
%----------------------------------------

figure;
t = (0:N-1)*dt;

plot(t, x(3,:), 'r', 'LineWidth', 1.5);
grid on;

xlabel('Time (s)');
ylabel('Height z (m)');
title('Droplet Height vs Time');


%----------------------------------------
    %the plot of spread%
%----------------------------------------

figure;

plot(x(1,:), x(2,:), 'k-', 'LineWidth', 1.5);
grid on;

xlabel('x (m)');
ylabel('y (m)');
title('Horizontal Spread of Droplet');
axis equal;








