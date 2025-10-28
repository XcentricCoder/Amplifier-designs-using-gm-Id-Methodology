%% ------------------------------------------------------------------------
%  Cascode Op-Amp Design using gm/ID Methodology
%  Works with nch_18.mat (BSIM data) and lookup.m (Murmann/Jespers format)
%  ------------------------------------------------------------------------
clear; clc; close all;

%% === Load transistor data ===
load('nch_18.mat');     % loads struct 'nch'

% Clean NaN and Inf from all relevant fields
fields_to_clean = {'GM','ID','GDS','CGG','CGD','CGS'};
for k = 1:length(fields_to_clean)
    f = fields_to_clean{k};
    if isfield(nch, f)
        tmp = nch.(f);
        tmp(~isfinite(tmp)) = 0;
        nch.(f) = tmp;
    end
end

%% === Input design targets ===
gm_ID_in  = 15;         % [S/A]  gm/ID of input transistor
gm_ID_cas = 10;         % [S/A]  gm/ID of cascode transistor
L_in      = 0.18e-6;    % [m]    input transistor channel length
L_cas     = 0.36e-6;    % [m]    cascode transistor channel length
I_total   = 100e-6;     % [A]    total tail current

% Clamp to valid gm/ID and L range
gm_ID_in  = max( min( gm_ID_in, 20 ), 2 );
gm_ID_cas = max( min( gm_ID_cas, 20 ), 2 );
L_in  = max(min(L_in, max(nch.L)), min(nch.L));
L_cas = max(min(L_cas, max(nch.L)), min(nch.L));

%% === Compute ratios using lookup() ===
% Each lookup computes Y/X vs gm/ID for given L

gm_gds_in  = lookup(nch, 'GM_GDS', 'GM_ID', gm_ID_in,  'L', L_in);
gm_cgg_in  = lookup(nch, 'GM_CGG', 'GM_ID', gm_ID_in,  'L', L_in);
gm_gds_cas = lookup(nch, 'GM_GDS', 'GM_ID', gm_ID_cas, 'L', L_cas);
gm_cgg_cas = lookup(nch, 'GM_CGG', 'GM_ID', gm_ID_cas, 'L', L_cas);

% Derived parameters
fT_in  = gm_cgg_in / (2*pi);   % fT = (gm/CGG)/(2π)
fT_cas = gm_cgg_cas / (2*pi);

%% === Compute small-signal parameters ===
% Split current equally between diff pair devices
ID_in  = I_total / 2;
gm_in  = gm_ID_in * ID_in;
gm_cas = gm_ID_cas * ID_in;

gds_in  = gm_in / gm_gds_in;
gds_cas = gm_cas / gm_gds_cas;

% Gain per branch (half circuit)
Av_in   = gm_in / gds_in;
Av_cas  = gm_cas / gds_cas;

% Cascode effective gain
Av0_branch = gm_in * (1/gds_in + 1/gds_cas);
Av0_total  = 2 * Av0_branch;   % differential output gain

% Unity-gain bandwidth (approx)
GBW = fT_in / (1 + gm_cas/gds_cas);

%% === Display results ===
fprintf('\n=== Cascode Op-Amp Design Results ===\n');
fprintf('gm/ID (input)      = %.2f V^-1\n', gm_ID_in);
fprintf('gm/ID (cascode)    = %.2f V^-1\n', gm_ID_cas);
fprintf('L_in               = %.2g m\n', L_in);
fprintf('L_cas              = %.2g m\n', L_cas);
fprintf('ID per branch      = %.2f µA\n', ID_in*1e6);
fprintf('fT_in              = %.2f MHz\n', fT_in/1e6);
fprintf('fT_cas             = %.2f MHz\n', fT_cas/1e6);
fprintf('|Av_in| (per branch)   = %.2f\n', Av_in);
fprintf('|Av_cas| (per branch)  = %.2f\n', Av_cas);
fprintf('|Av0_total| (diff)     = %.2f\n', Av0_total);
fprintf('GBW (approx)           = %.2f MHz\n', GBW/1e6);
fprintf('------------------------------------\n');

%% === Optional: plot gm/gds and gm/cgg vs gm/ID ===
gm_ID_range = linspace(2,20,50);
gm_gds_plot = arrayfun(@(x) lookup(nch,'GM_GDS','GM_ID',x,'L',L_in), gm_ID_range);
gm_cgg_plot = arrayfun(@(x) lookup(nch,'GM_CGG','GM_ID',x,'L',L_in), gm_ID_range);

figure;
subplot(1,2,1);
plot(gm_ID_range, gm_gds_plot, 'LineWidth', 1.5);
xlabel('gm/ID [S/A]'); ylabel('gm/gds');
title('Intrinsic Gain vs gm/ID'); grid on;

subplot(1,2,2);
plot(gm_ID_range, gm_cgg_plot/(2*pi), 'LineWidth', 1.5);
xlabel('gm/ID [S/A]'); ylabel('f_T [Hz]');
title('Transition Frequency vs gm/ID'); grid on;
%% === Plot Gain vs Frequency (Bode Magnitude) ===
% --------------------------------------------------
f = logspace(3, 9, 500);             % frequency range: 1 kHz to 1 GHz
fp = GBW / Av0_total;                % dominant pole frequency
Av_f = Av0_total ./ sqrt(1 + (f/fp).^2);  % single-pole gain response

figure;
semilogx(f, 20*log10(Av_f), 'k', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('|A_v| (dB)');
title('Open-Loop Gain vs Frequency (Cascode Op-Amp)');
grid on;
text(fp, 20*log10(Av0_total)-3, sprintf('f_p = %.2f MHz', fp/1e6), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
text(GBW, 0, sprintf('GBW = %.2f MHz', GBW/1e6), ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
