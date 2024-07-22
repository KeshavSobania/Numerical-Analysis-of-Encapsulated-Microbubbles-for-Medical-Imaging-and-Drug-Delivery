% Given parameters
R0 = 1.25e-6; % Initial bubble radius normalized to 1
patm = 101325; % Atmospheric pressure in Pa
kA = 2.05e-9; % Coefficient of diffusivity of air in water in m^2/s
hA = 2.857e-5; % Permeability of air through encapsulation in m/s
c0 = 0.025; % Surface tension in N/m
LA = 1.71e-2; % Ostwald coefficient of OFP for LA
f=1;
% Values of es to consider
es_values = [c0, 1.05*c0, 1.2*c0, 1.5*c0, 2*c0, 4*c0, 20*c0];

% Time parameters
t_start = 1e-3; % Initial time in seconds
t_end = 1e3; % End time in seconds
num_steps = 1000000; % Number of time steps
t = logspace(log10(t_start), log10(t_end), num_steps);

% Plotting the dissolution dynamics for different es values
figure;
hold on;
for i = 1:length(es_values)
    es = es_values(i);
    ag = kA / (hA * R0);
    e = 2 * es / (patm * R0);

    % Calculate dimensionless time parameter tau
    tao = (t * kA) / R0^2;

    % Initial condition
    r = ones(1, num_steps); % r(t=0) = 1

    % Numerical integration using provided equations
    for j = 2:num_steps
        y = 2 * c0 / (patm * R0);

        if y + e * (r(j-1)^2 - 1) < 0
            dr_dt = (kA/R0^2)*(-3 * LA / (ag + r(j-1))) * ((r(j-1) * (1 - f)) / (3 * r(j-1)));
        else
            dr_dt = (kA/R0^2)*(-3 * LA / (ag + r(j-1))) * ((r(j-1) * (1 - f) + y + e * (r(j-1)^2 - 1)) / (3 * r(j-1) + 2 * y + 2 * e * (2 * r(j-1)^2 - 1)));
        end

        % Update r using Runge-Kutta method
        h = t(j) - t(j - 1);
        k1 = h * dr_dt;
        k2 = h * (-3 * LA / (ag + (r(j-1) + k1/2))) * ((r(j-1) * (1 - 1)) / (3 * r(j-1)));
        k3 = h * (-3 * LA / (ag + (r(j-1) + k2/2))) * ((r(j-1) * (1 - 1)) / (3 * r(j-1)));
        k4 = h * (-3 * LA / (ag + (r(j-1) + k3))) * ((r(j-1) * (1 - 1)) / (3 * r(j-1)));

        r(j) = r(j - 1) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

    % Plotting the results
    semilogx(t, r, 'DisplayName', ['es = ' num2str(es/c0) 'c0']);
end

xlabel('Time (s)');
ylabel('r (R/R0)');
title('Dissolution Dynamics of Encapsulated Air Bubble for f=1 saturated air medium');
legend('show');
grid on;
hold off;
