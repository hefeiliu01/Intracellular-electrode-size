clear all

% Assumption: cells are tightly arranged in a hexagonal way.
% Electode position uniformly distributed in the lattice

% Set electrode size paraters
rcell = 10e-6; % unit: m
min_ratio = 1e-4;
max_ratio = 3;

% Initialize the electrode radius vector, electrode impedance vector
num_size = 101;
relectrode = logspace(log10(min_ratio * rcell), log10(max_ratio * rcell), num_size);
electrode_area = pi * relectrode.^2;

height = 2e-6; % pillar height 2 um
electrode_surface_area = 2 * pi * relectrode * height;

% Set circiut parameters
Rs = 0.1e9; %unit: ohm
Rjm = 400e6; %unit: ohm
Zin = 50e6; %unit: ohm

% Set Ze0, using SH5, Ze = 300 kOhm at 5 kHz, or 1.5 MOhm at 1 kHz
r0 = 5e-6; %unit: m
rin0 = 4e-6; %unit: m
hole_depth = 1.4e-6; %unit: m
Area0 = pi * r0^2 + 2 * pi * rin0 * hole_depth;
Z0 = 1.5e6; %unit: ohm, at 1 kHz
Y0 = 1 / Z0;
Ye0 = Y0 / Area0;
Ze0 = Z0 * Area0 / 1; % 1 for normal Ze0, 10 for reduced Ze0

Ze0 = 50e-6; % set from 5e-6 to 500e-6

% Define the function for parallel impedance calculations
function Z = impedance_para(Z1, Z2)
    if (Z1 == Inf)
        Z = Z2;
        else if (Z1 == Inf)
            Z = Z1;
            else
            Z = (Z1 * Z2) / (Z1 + Z2);
        end
    end
end

% Set the final output, mean alpha
mean_alpha = zeros(1,num_size);




% Start Monte-Carlo sampling the 2d space, calculate the overlap area
% Define the coordinates for the sampling grid
num_step = 101;

x = linspace(0, rcell, num_step);
y = linspace(0, rcell/sqrt(3), num_step);

for k = 1:num_size

    % Initialize variables for storing results
    num_effect_point = 0;
    sum_alpha = 0;

    % Initialize the R and r of two circles for different relectrode and rcell
    R = max(rcell,relectrode(k));
    r = min(rcell,relectrode(k));
    D1 = R - r;
    D2 = R + r;
    area_small_circle = pi * r.^2;

    % set the electrode parameters

    % simply consider the sidewall height = 4 * re
    % the surface area fo sidewall = 8 * pi * re^2

    Ze = Ze0 / electrode_surface_area(k); % nanopillar impedance

for i = 1:num_step
    for j = 1:num_step
        if j <= i
            num_effect_point = num_effect_point + 1;

        %calculate phi
            d = sqrt(x(i)^2 + y(j)^2);
            if d <= D1
                overlap_area = area_small_circle;
            end
            
            if (d > D1) & (d <= D2)
                theta1 = acos((d^2 + R^2 - r^2) / (2 * d * R));
                theta2  = acos((d^2 + r^2 - R^2) / (2 * d * r));

                overlap_area = R^2 * theta1 + r^2 * theta2 ...
                            - 0.5 * sqrt( ...
                            (-d + R + r) * ...
                            ( d + R - r) * ...
                            ( d - R + r) * ...
                            ( d + R + r) );
            end
            if d > D2
                overlap_area = 0;
            end

            phi = overlap_area / electrode_area(k);

        %calculate alpha
            Zee = (1 / phi) * Ze;
            Zep = (1 / (1 - phi)) * Ze;

            Zine = Zee + impedance_para(Zep,Zin);
            
            tran1 = impedance_para(Zine,Rs) / (impedance_para(Zine,Rs) + Rjm);
            tran2 = impedance_para(Zep,Zin) / (Zee + impedance_para(Zep,Zin));
            alpha = tran1 * tran2;

            sum_alpha = sum_alpha + alpha; % Accumulate the product of x values
        end
    end
end

mean_alpha(k) = sum_alpha / num_effect_point;

end

relectrode_tran = relectrode';
delectrode_tran = 2 * relectrode_tran;
mean_alpha_tran = mean_alpha';

% output meanalpha vs. size
% meanAlpha = table(2 * relectrode', mean_alpha', 'VariableNames', {'ElectrodeDiameter (um)', 'MeanAlpha'});
% writetable(meanAlpha, 'Disk MeanAlpha vs ElectrodeDiameter.csv');

% plot meanalpha vs. size
figure(1)
plot(2 * relectrode, mean_alpha)
xlabel('electrode diameter (m), cell diameter = 20 um')
ylabel('mean alpha')
set(gca, 'XScale', 'log')
%title('Ze ~ 1/A')
title('Ze0=5 MOhm*um^2, Rs = 100 MOhm, varying Zjm, pillar height 2 um')
legend('400 MOhm','100 MOhm','25 MOhm','1 MOhm','1 kOhm')
hold on;

