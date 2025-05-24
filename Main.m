clear
close all
clc

path = pwd;        
cd Data
%% Enter Structural Information
Samplename = 'H2O';                      % Specify the sample for refractive index calculation
Reference = load('H2O_Reference.txt');   % Load the experimental time trace of the reference
Sample = load('H2O_Sample.txt');         % Load the experimental time trace of the sample
PopSize = 100;                           % Population size: higher gives better accuracy but increases run time
Maxit = 1000;                            % Max generations: more iterations improve results but take longer
data = jsondecode(fileread([Samplename,'.json']))                        
nr = data.("settings").calibration_index;                  
minTHz = data.("settings").minTHz;                        
maxTHz = data.("settings").maxTHz;                                                
d = [data.Sample.d_nm];                                    
nk = [data.Sample.n];
t_smpl0 = d(find(isnan(nk)));                              
                          
 
%% Calculate experimental Transfer Function and Fourier Transform
cd(path)                              
t_file_ref = Reference(:,1) * 1e-12;        
E_file_ref = Reference(:,2);               
t_file_sig = Sample(:,1) * 1e-12;       
E_file_sig = Sample(:,2);               
Time2Freq(t_file_ref,E_file_ref, t_file_sig, E_file_sig, minTHz, maxTHz)

%% Analytical extraction of n and k
clearvars -except ii dataname t_smpl0 nr d nk PopSize Maxit Samplename % Preserve key vars
load("Test.mat")
c = 299792458;                              % Speed of light (m/s)
L = numel(lambda0);                         % Number of wavelengths
nn0 = 1 + c * dtpeaks / (t_smpl0 * 1e-9);   % closEffective refractive index from time delay
neff = nk; neff(isnan(neff)) = nn0;         % Update refractive indices (replace NaNs with estimated nn0)
dlimit = (c * dT ./ (2 * neff)) * 1e9;      % Threshold thickness for each layer
PH = 0:5;                                   % First 6 unwrapped phases (the correct one is 0)
delta_phi_values = 2*pi*floor((PH+1)/2) .* (-1).^PH;

% Loop for calculating analytical refractive index from different unwrapped phases
for i = 1:length(delta_phi_values)
    delta_phi2 = delta_phi + delta_phi_values(i);
    n_anlt = nr + c * delta_phi2 ./ (2 * pi * f.' * t_smpl0 * 1e-9);
    constnt = (4 .* n_anlt .* nr) ./ ((abs(EsovEr) .* (n_anlt + nr).^2));
    k_anlt = (c ./ (2 * pi * f.' * t_smpl0 * 1e-9)) .* log(constnt);    
    % Store results
    n_anltic(i,:) = n_anlt.';
    k_anltic(i,:) = k_anlt.';
end
% plot(f, n_anltic(1,:));
% hold on
% plot(f, k_anltic(1,:));
save('Test.mat')  % Save updated analytical results
%% Optimization Process
%Define Optimization Bounds from Analytical Estimates
nhalf = [(n_anltic(1,:) + n_anltic(2,:)) / 2;
         (n_anltic(1,:) + n_anltic(3,:)) / 2];          
khalf = -[(k_anltic(1,:) + k_anltic(2,:)) / 2;
          (k_anltic(1,:) + k_anltic(3,:)) / 2];        
lb = [min(real(nhalf)), -2*max(abs(k_anltic(:))).*ones(1,L), t_smpl0];  % Lower bounds
ub = [max(real(nhalf)), zeros(1,L), t_smpl0];  % Upper bounds
% figure
% plot(f, lb(L+1:2*L));
% hold on
% plot(f, ub(L+1:2*L));
% Plotting options (not used here, but defined for clarity)
plot_opts = {'LineStyle', ':', 'Marker', 'o', 'LineWidth', 1.6};
axis_opts = {'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5};
fig_opts = {'Units', 'Inches', 'Position', [1, 1, 6, 4]};

% Optional: plot bounds for n and k
% figure; plot(f, lb(1:L), plot_opts{:}); hold on; plot(f, ub(1:L), plot_opts{:});
% figure; plot(f, lb(L+1:2*L), plot_opts{:}); hold on; plot(f, ub(L+1:2*L), plot_opts{:});

%% Run Genetic Algorithm
nvars = numel(lb);                         % Number of variables
warning('off','all')                       % Suppress warnings
fun = @TM_DBR_test1;                       % Objective function
initialPop = [n_anltic(1,:), -k_anltic(1,:), t_smpl0];  % Initial guess
options = optimoptions('ga', ...
    'PopulationSize', PopSize, ...
    'InitialPopulationMatrix', initialPop, ...
    'MaxGenerations', Maxit, ...
    'UseParallel', true, ...
    'Display', 'iter', ...
    'PlotFcn', @gaplotbestf);              % GA settings
[d0, fval, exitflag, output] = ga(fun, nvars, [], [], [], [], lb, ub, [], options);
cd Results
save([Samplename, '_Results']);                 % Save result

%% Plot the results
f=(c./lambda0)*1e-3;
% Plot refractive index (n)
figure
subplot(1, 2, 1);
plot(f, d0(1:L), plot_opts{:}); %hold on;
%plot(f, n_anltic(1,:), plot_opts{:});
ylabel('Refractive index, n', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Frequency (THz)', 'FontSize', 12, 'FontWeight', 'bold')
%legend({'MTMM', 'Analyt._{\Delta\phi}'},'Location', 'northoutside', 'FontSize', 12, 'Orientation', 'horizontal')
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'Box', 'on','FontWeight', 'bold')
set(gca, 'LineWidth', 1.5)
%set(gcf, 'Position', [100, 100, 600, 400])
% Plot extinction coefficient (k)
subplot(1, 2, 2);
plot(f, -d0(L+1:2*L), plot_opts{:}); %hold on;
%plot(f, k_anltic(1,:), plot_opts{:});
ylabel('Extinction coefficient, k', axis_opts{:});
xlabel('Frequency (THz)', 'FontSize', 12, 'FontWeight', 'bold')
%legend({'MTMM', 'Analyt._{\Delta\phi}'},'Location', 'northoutside', 'FontSize', 12, 'Orientation', 'horizontal')
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'Box', 'on','FontWeight', 'bold')
set(gca, 'LineWidth', 1.5)
set(gcf, 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.3]);
saveas(gcf, [Samplename, '_Results.fig']);


function objctv = TM_DBR_test1(d0)
%% Load Optimization Results and Data
load("Test.mat");                          % Load frequency-dependent variables 
L = numel(lambda0);                        % Number of wavelengths
[ns(1:L), ksmp(1:L)] = deal(d0(1:L), d0(L+1:2*L));  % Extract real (n) and imaginary (k) parts
ns = ns + 1i * ksmp;                      % Combine to form complex refractive index
t_smpl = d0(2*L + 1);                     % Extract sample thickness
theta0 = 0;                               % Normal incidence
idx = find(isnan(d));                     % Find sample layer (NaN)
d(idx) = t_smpl;                          % Insert sample thickness
flag = 0;                                 % Reference calculation
t_cs2 = MTMM(d, lambda0, theta0, nr, ns, flag, dlimit, nk);
flag = 1;                                 % Sample calculation
t_cs3 = MTMM(d, lambda0, theta0, nr, ns, flag, dlimit, nk);
deviations = abs(EsovEr - (t_cs3 ./ t_cs2));  % Deviation from experimental field ratio
objctv = sum(deviations);                    % Sum of deviations (objective value)
end
