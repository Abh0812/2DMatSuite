%% Generic Ellipsometry Analysis for a 2D Material on SiO₂/Si
% This script performs a point-by-point inversion to extract the 
% refractive index (n) and extinction coefficient (k) of a 2D material 
% deposited on a SiO₂/Si substrate.
%
% INPUTS:
%   - SiO₂ data file (e.g., "SiO2.txt"): two columns [Energy (eV), n_SiO2]
%       (Assumes k_SiO2 is negligible and set to 0)
%   - Si data files:
%       * Si n file (e.g., "Si_n.txt"): two columns [Energy (eV), n_Si]
%       * Si k file (e.g., "Si_k.txt"): two columns [Energy (eV), k_Si]
%   - Experimental data file (e.g., "exp_data.txt"): three columns 
%       [Wavelength (nm), Δ (radians), Ψ (radians)]
%
%   - Film (2D material) thickness (nm) and SiO₂ thickness (nm)
%   - Angle of incidence (degrees)
%
% The script automatically sets the wavelength range for the analysis 
% based on the minimum and maximum wavelengths in the experimental file.

clear; clc;

%% User Inputs: File Names & Experimental Parameters
sio2_filename = input('Enter the SiO₂ data filename (e.g., "SiO2.txt"): ', 's');
si_n_filename = input('Enter the Si n data filename (e.g., "Si_n.txt"): ', 's');
si_k_filename = input('Enter the Si k data filename (e.g., "Si_k.txt"): ', 's');
exp_filename  = input('Enter the experimental data filename (e.g., "exp_data.txt"): ', 's');

d_film = input('Enter the film thickness (nm): ');
d_SiO2 = input('Enter the SiO₂ thickness (nm): ');
theta_inc_deg = input('Enter the angle of incidence (degrees): ');
theta_inc = deg2rad(theta_inc_deg);

%% Load Known Refractive Index Data for SiO₂
% File format: two columns [Energy (eV), n_SiO2]
sio2_data = load(sio2_filename);
energy_sio2 = sio2_data(:,1);
n_sio2 = sio2_data(:,2);
k_sio2 = zeros(size(n_sio2));  % k_SiO2 assumed negligible
wl_sio2 = 1240 ./ energy_sio2;  % Convert energy (eV) to wavelength (nm)

% (Optional) If desired, restrict to a specific wavelength range.
mask_sio2 = (wl_sio2 >= min(wl_sio2)) & (wl_sio2 <= max(wl_sio2));
wl_sio2 = wl_sio2(mask_sio2);
n_sio2 = n_sio2(mask_sio2);
k_sio2 = k_sio2(mask_sio2);
[wl_sio2, sortIdx] = sort(wl_sio2);
n_sio2 = n_sio2(sortIdx);
k_sio2 = k_sio2(sortIdx);
n_SiO2_interp = @(wl) interp1(wl_sio2, n_sio2, wl, 'linear', 'extrap');
k_SiO2_interp = @(wl) interp1(wl_sio2, k_sio2, wl, 'linear', 'extrap');

%% Load Known Refractive Index Data for Si
% Si n file: [Energy (eV), n_Si]
si_n_data = load(si_n_filename);
energy_si_n = si_n_data(:,1);
n_si = si_n_data(:,2);
% Si k file: [Energy (eV), k_Si]
si_k_data = load(si_k_filename);
energy_si_k = si_k_data(:,1);
k_si = si_k_data(:,2);
wl_si_n = 1240 ./ energy_si_n;
wl_si_k = 1240 ./ energy_si_k;
mask_si_n = (wl_si_n >= min(wl_si_n)) & (wl_si_n <= max(wl_si_n));
mask_si_k = (wl_si_k >= min(wl_si_k)) & (wl_si_k <= max(wl_si_k));
wl_si_n = wl_si_n(mask_si_n);
n_si = n_si(mask_si_n);
wl_si_k = wl_si_k(mask_si_k);
k_si = k_si(mask_si_k);
[wl_si_n, idx_n] = sort(wl_si_n);
n_si = n_si(idx_n);
[wl_si_k, idx_k] = sort(wl_si_k);
k_si = k_si(idx_k);
n_Si_interp = @(wl) interp1(wl_si_n, n_si, wl, 'linear', 'extrap');
k_Si_interp = @(wl) interp1(wl_si_k, k_si, wl, 'linear', 'extrap');

%% Load Experimental Ellipsometry Data
% File format: three columns [Wavelength (nm), Δ (rad), Ψ (rad)]
exp_data = load(exp_filename);
wl_exp = exp_data(:,1);
delta_exp = exp_data(:,2);
psi_exp = exp_data(:,3);

% Automatically determine wavelength range from experimental data
wl_min = min(wl_exp);
wl_max = max(wl_exp);
fprintf('Experimental wavelength range: %.1f nm to %.1f nm\n', wl_min, wl_max);

%% Inversion to Extract Film n and k (Point-by-Point)
% Define bounds for the film parameters: [n, k]
lb = [0.1, 0.0];
ub = [5.0, 5.0];

opts = optimoptions('fmincon', 'Display', 'off');
n_extracted = zeros(length(wl_exp), 1);
k_extracted = zeros(length(wl_exp), 1);

% Loop over each experimental wavelength
for i = 1:length(wl_exp)
    wl = wl_exp(i);
    psi_meas = psi_exp(i);
    delta_meas = delta_exp(i);
    
    % Objective function: squared error between modeled and measured Ψ and Δ
    objFun = @(x) objective_fun(x, wl, theta_inc, psi_meas, delta_meas, ...
                                 d_film, d_SiO2, n_SiO2_interp, k_SiO2_interp, n_Si_interp, k_Si_interp);
    x0 = [1.5, 0.1];  % Initial guess for [n, k]
    
    [x_opt, ~] = fmincon(objFun, x0, [], [], [], [], lb, ub, [], opts);
    n_extracted(i) = x_opt(1);
    k_extracted(i) = x_opt(2);
end

%% Plot the Results
figure;
subplot(1,2,1);
plot(wl_exp, n_extracted, 'bo-','LineWidth',1.5);
xlabel('Wavelength (nm)');
ylabel('Refractive Index, n');
title('Extracted n vs. Wavelength');
grid on;

subplot(1,2,2);
plot(wl_exp, k_extracted, 'ro-','LineWidth',1.5);
xlabel('Wavelength (nm)');
ylabel('Extinction Coefficient, k');
title('Extracted k vs. Wavelength');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = objective_fun(x, wl, theta_inc, psi_meas, delta_meas, ...
                             d_film, d_SiO2, n_SiO2_interp, k_SiO2_interp, n_Si_interp, k_Si_interp)
    % x contains [n_film, k_film] for the 2D material.
    n_film = x(1);
    k_film = x(2);
    [psi_model, delta_model] = calc_multilayer_psi_delta(n_film, k_film, wl, theta_inc, ...
                                                          d_film, d_SiO2, n_SiO2_interp, k_SiO2_interp, n_Si_interp, k_Si_interp);
    err = (psi_model - psi_meas)^2 + (delta_model - delta_meas)^2;
end

function [psi, delta] = calc_multilayer_psi_delta(n_film, k_film, wl, theta_inc, d_film, d_SiO2, n_SiO2_interp, k_SiO2_interp, n_Si_interp, k_Si_interp)
    % Calculate the ellipsometric parameters Ψ and Δ for the stack:
    % Air / 2D Material (film) / SiO₂ / Si
    
    N0 = 1.0;  % Air
    N_film = n_film - 1i*k_film;  % Film (2D material)
    
    % SiO₂ layer properties (from interpolation)
    n_SiO2 = n_SiO2_interp(wl);
    k_SiO2 = k_SiO2_interp(wl);
    N_SiO2 = n_SiO2 - 1i*k_SiO2;
    
    % Si substrate properties (from interpolation)
    n_Si = n_Si_interp(wl);
    k_Si = k_Si_interp(wl);
    N_Si = n_Si - 1i*k_Si;
    
    % Compute refraction angles using Snell's law (with real parts)
    theta0 = theta_inc;
    theta_film = asin((N0 * sin(theta0)) / real(N_film));
    theta_SiO2 = asin((N0 * sin(theta0)) / real(N_SiO2));
    theta_Si = asin((N0 * sin(theta0)) / real(N_Si));
    
    % Vacuum wavenumber (nm^-1)
    k0 = 2*pi / wl;
    
    % Phase thicknesses
    delta_film = k0 * N_film * d_film * cos(theta_film);
    delta_SiO2 = k0 * N_SiO2 * d_SiO2 * cos(theta_SiO2);
    
    % Fresnel reflection coefficients for s and p polarizations
    r_s = @(n1, n2, th1, th2) (n1*cos(th1) - n2*cos(th2)) / (n1*cos(th1) + n2*cos(th2));
    r_p = @(n1, n2, th1, th2) (n2*cos(th1) - n1*cos(th2)) / (n2*cos(th1) + n1*cos(th2));
    
    % Calculate reflection coefficients at interfaces
    r01_s = r_s(N0, N_film, theta0, theta_film);
    r01_p = r_p(N0, N_film, theta0, theta_film);
    r12_s = r_s(N_film, N_SiO2, theta_film, theta_SiO2);
    r12_p = r_p(N_film, N_SiO2, theta_film, theta_SiO2);
    r23_s = r_s(N_SiO2, N_Si, theta_SiO2, theta_Si);
    r23_p = r_p(N_SiO2, N_Si, theta_SiO2, theta_Si);
    
    % Propagation factors for film and SiO₂ layers
    P_film = exp(1i * delta_film);
    P_SiO2 = exp(1i * delta_SiO2);
    
    % Combine layers: first combine SiO₂ and Si into an effective reflection
    r_eff_12_s = (r12_s + r23_s * P_SiO2^2) / (1 + r12_s * r23_s * P_SiO2^2);
    r_eff_12_p = (r12_p + r23_p * P_SiO2^2) / (1 + r12_p * r23_p * P_SiO2^2);
    
    % Combine with the top interface (Air/Film)
    r_total_s = (r01_s + r_eff_12_s * P_film^2) / (1 + r01_s * r_eff_12_s * P_film^2);
    r_total_p = (r01_p + r_eff_12_p * P_film^2) / (1 + r01_p * r_eff_12_p * P_film^2);
    
    % Ellipsometric parameters from the ratio ρ = r_total_p / r_total_s
    rho = r_total_p / r_total_s;
    psi = atan(abs(rho));
    delta = angle(rho);
end
