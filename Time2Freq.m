function [] = Time2Freq(t_ref,E_ref,t_sam,E_sam,minTHz,maxTHz)
t_ref = t_ref - min(t_ref);
t_sam = t_sam - min(t_sam);

figure
plot(t_ref,E_ref, 'lineWidth',1.5)
hold on
plot(t_sam,E_sam, 'lineWidth',1.5)
xlabel('Time(sec)','FontSize',16)
ylabel('Electric field intensity ({a.u.})','linewidth',16,'FontSize',16);
legend('E_{Reference}', 'E_{Sample}')

t_min = min([min(t_ref), min(t_sam)]);
t_max = max([max(t_ref), max(t_sam)]);

% Define a new common time grid
num_points = max(length(t_ref), length(t_sam));
time = linspace(t_min, t_max, num_points)';

% Interpolate signals to common time grid (optional, can uncomment if needed)
% E_ref = interp1(t_ref, E_ref, time, 'linear', 0);
% E_sam = interp1(t_sam, E_sam, time, 'linear', 0);

% === Frequency bounds in Hz ===
f_min = minTHz * 1e12;
f_max = maxTHz * 1e12;

% === Step 1: Locate pulse maxima ===
[~, idx_ref] = max((E_ref));
[~, idx_sam] = max((E_sam));
t0r = time(idx_ref);
t0s = time(idx_sam);
dtpeaks = t0s - t0r;
dtmin = dtpeaks;
dT = t_max - t_min;

% === Step 2: Compute one-sided FFT with zero-padding ===
N = length(time);
dt = time(2) - time(1);
Fs = 1/dt;

pad_factor = 1;
N_pad = pad_factor * N;

% Zero-padding
E_ref_padded = [E_ref; zeros(N_pad - N, 1)];
E_sam_padded = [E_sam; zeros(N_pad - N, 1)];

% Frequency vector
f_full = Fs * (0:floor(N_pad/2)) / N_pad;
omega_full = 2 * pi * f_full;

% FFTs
E_ref_fft_full = fft(E_ref_padded);
E_sam_fft_full = fft(E_sam_padded);

% One-sided spectra
E_ref = E_ref_fft_full(1:length(f_full));
E_sam = E_sam_fft_full(1:length(f_full));

% === Step 3: Filter frequency range ===
freq_mask = (f_full >= f_min) & (f_full <= f_max);
f = f_full(freq_mask);
omega = omega_full(freq_mask);
E_ref = E_ref(freq_mask);
E_sam = E_sam(freq_mask);

% === Step 4: Compute reduced phase ===
phi0_ref = omega * t0r;
phi0_sam = omega * t0s;

phi_red_ref = angle(E_ref .* exp(-1i * phi0_ref.'));
phi_red_sam = angle(E_sam .* exp(-1i * phi0_sam.'));

% === Step 5: Unwrap phase difference ===
delta_phi_star_0 = unwrap(phi_red_sam - phi_red_ref);

% === Step 6: Linear fit and offset removal ===
center_fraction = 0.5;
N_center = round(length(f) * center_fraction);
start_idx = round((length(f) - N_center)/2);
center_idx = start_idx : start_idx + N_center - 1;

omega_center = omega(center_idx);
delta_phi_center = delta_phi_star_0(center_idx);

p = polyfit(omega_center, delta_phi_center, 1);
B = p(2);
delta_phi_0 = delta_phi_star_0 - 2*pi * round(B / (2*pi));

% === Step 7: Final corrected phase ===
phi_offset = 0;
delta_phi = -1 * (delta_phi_0 - phi0_ref.' + phi0_sam.' + phi_offset);

% === Plotting ===
figure;
plot(f, log10(abs(E_ref)), 'lineWidth', 1.5);
hold on
plot(f, log10(abs(E_sam)), 'lineWidth', 1.5);
title('One-sided Fourier Transform');
xlabel('Frequency(Hz)', 'FontSize', 14)
ylabel('Electric field intensity ({a.u.})', 'FontSize', 14)
legend('E_{Reference}', 'E_{Sample}')
ax = gca;
ax.FontSize = 12;

figure;
plot(f * 1e-12, delta_phi, 'k', 'LineWidth', 1.5);
xlabel('Frequency (THz)');
ylabel('Phase Difference (rad)');
title('Corrected Phase Difference');
grid on;

% === Save Data ===
c = 299792458; % Speed of light (m/s)
f = flip(f);
EsovEr = flip(E_sam ./ E_ref);
E_sam = flip(E_sam);
E_ref = flip(E_ref);
delta_phi = flip(delta_phi);
lambda0 = (c ./ f) * 1e9;

save('Test.mat', 'EsovEr', 'f', 'E_sam', 'dtmin', 'dT', 'E_ref', 'lambda0', 'dtpeaks', 't0s', 't0r', 'delta_phi');

end
