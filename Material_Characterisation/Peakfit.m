%% Peak Fitting Script for Raman/PL Spectrum using Gaussian or Lorentzian Models
% Clear workspace and command window
clear; clc;

%% Step 1: Load the Data
% Assumes the data file (e.g., 'spectrum.txt') has two columns:
%   column 1: wavenumber or wavelength
%   column 2: intensity
dataFile = input('Enter the data filename (e.g., spectrum.txt): ', 's');
data = load(dataFile);
x = data(:,1); % wavenumber/wavelength
y = data(:,2); % intensity

%% Step 2: Get User Inputs
% Ask for model type
modelType = input('Enter model type ("gaussian" or "lorentzian"): ', 's');

% Ask for number of peaks to fit
numPeaks = input('Enter the number of peaks to fit: ');

% Instead of prompting one-by-one, ask for comma-separated arrays.
fprintf('\nEnter the lower and upper bounds as MATLAB arrays (e.g., [val1, val2, ...]).\n');

% For the peak centers:
centerLower = input('Enter lower bounds for peak centers: ');
centerUpper = input('Enter upper bounds for peak centers: ');

% For the FWHM:
fwhmLower   = input('Enter lower bounds for peak FWHM: ');
fwhmUpper   = input('Enter upper bounds for peak FWHM: ');

% Check that the entered arrays have the correct number of elements.
if length(centerLower) ~= numPeaks || length(centerUpper) ~= numPeaks || ...
   length(fwhmLower) ~= numPeaks || length(fwhmUpper) ~= numPeaks
    error('Number of bounds entered does not match the number of peaks specified.');
end

%% Step 3: Set Up Initial Guesses and Bounds for lsqcurvefit
% Each peak is described by three parameters: [Amplitude, Center, FWHM]
p0 = []; lb = []; ub = [];
for i = 1:numPeaks
    % Estimate amplitude from the maximum intensity near the provided center bounds
    idx = find(x >= centerLower(i) & x <= centerUpper(i));
    if isempty(idx)
        ampGuess = max(y);
    else
        ampGuess = max(y(idx));
    end
    centerGuess = (centerLower(i) + centerUpper(i)) / 2;
    fwhmGuess   = (fwhmLower(i) + fwhmUpper(i)) / 2;
    
    p0 = [p0; ampGuess; centerGuess; fwhmGuess];
    lb = [lb; 0; centerLower(i); fwhmLower(i)];
    ub = [ub; Inf; centerUpper(i); fwhmUpper(i)];
end

%% Step 4: Define the Composite Model Function
% The composite model sums the contributions from all peaks.
% For a Gaussian, the function is:
%    G(x) = A * exp( -4*log(2)*((x-center)/FWHM).^2 )
% For a Lorentzian, the function is:
%    L(x) = A * ((0.5*FWHM)^2 ./ ((x-center).^2 + (0.5*FWHM)^2) )
modelFun = @(p, xdata) compositePeakModel(p, xdata, numPeaks, modelType);

%% Step 5: Perform the Curve Fit Using lsqcurvefit
options = optimset('Display','iter');
[p_fit, resnorm, residual, exitflag, output] = lsqcurvefit(modelFun, p0, x, y, lb, ub, options);

% Evaluate the fitted model
y_fit = modelFun(p_fit, x);

%% Step 6: Extract and Display Fitted Parameters
fprintf('\nFitted Peak Parameters:\n');
for i = 1:numPeaks
    A      = p_fit(3*i - 2);
    center = p_fit(3*i - 1);
    fwhm   = p_fit(3*i);
    fprintf('Peak %d: Amplitude = %.3f, Center = %.3f, FWHM = %.3f\n', i, A, center, fwhm);
end

%% Step 7: Plot Experimental Data and Fitted Curve
figure;
plot(x, y, 'b.', 'DisplayName', 'Experimental Data');
hold on;
plot(x, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Model');
xlabel('Wavenumber / Wavelength');
ylabel('Intensity');
legend('show');
title(sprintf('Peak Fitting using %s Model', modelType));
grid on;

%% Local Function: Composite Peak Model
function ysum = compositePeakModel(p, xdata, numPeaks, modelType)
    % Initialize the sum of peaks
    ysum = zeros(size(xdata));
    for i = 1:numPeaks
        % Extract parameters for the i-th peak
        A      = p(3*i - 2);
        center = p(3*i - 1);
        FWHM   = p(3*i);
        
        switch lower(modelType)
            case 'gaussian'
                % Gaussian peak: A * exp(-4*log(2)*((x-center)/FWHM).^2)
                ysum = ysum + A * exp(-4*log(2)*((xdata - center)/FWHM).^2);
            case 'lorentzian'
                % Lorentzian peak: A * ((0.5*FWHM)^2 ./ ((x-center).^2 + (0.5*FWHM)^2))
                ysum = ysum + A * ((0.5*FWHM)^2 ./ ((xdata - center).^2 + (0.5*FWHM)^2));
            otherwise
                error('Unknown model type: %s. Choose "gaussian" or "lorentzian".', modelType);
        end
    end
end
