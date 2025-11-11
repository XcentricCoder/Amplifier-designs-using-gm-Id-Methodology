% --------------------------
load nch_18.mat;
load pch_18.mat;
% Parameters (edit as needed)
fu = 1e9;            % target unity-gain frequency (Hz)
CL = 1e-12;          % load capacitance (F)
FO  = 10;            % fan-out (=> fT = fu*FO)
fT  = fu*FO;         % transistor fT requirement (Hz)
A_boost = 20;        % assumed gain of the local booster (linear). try 10, 20, 50...
gm_ID_load = 10;     % (gm/ID) for p-load (S/A), typical compromise as in book.
% L sweeps
L1 = 0.18:0.001:0.4;             % candidate gate lengths for M1 (um)
L_cas = [ 0.18, 0.3, 0.5, 1.0]; % sweep for cascode gate length (um)

% Convert to base units for calculation (if the lookup tables expect them in meters, 
% or keep as um if lookup tables are designed for um, but print in um)
L1 = L1 * 1e-6; % Converting to meters (assuming lookup tables use base units)
L_cas = L_cas * 1e-6; % Converting to meters

% 1) find gm/ID for M1 required to meet fT at each L1 (use GM_CGG->fT lookup)
gm_ID1_vec = lookup(nch,'GM_ID','GM_CGG',2*pi*fT,'L',L1);   % vector vs L1
gds_ID1_vec = diag( lookup(nch,'GDS_ID','GM_ID',gm_ID1_vec,'L',L1) );

% 2) for each cascode length compute its gds/ID (for chosen gm/ID_cas)
gm_ID_cas = 10;   % pick an initial gm/ID for the cascode device (you may sweep this)
gds_ID_cas = zeros(length(L1), length(L_cas));
for k = 1:length(L_cas)
    % This lookup returns a single value, replicated to match the size of L1 for vectorized division later.
    gds_ID_cas(:,k) = lookup(nch,'GDS_ID','GM_ID',gm_ID_cas,'L',L_cas(k)); % n-type for cascode
end

% 3) incorporate booster effect: reduce cascode conductance by (1 + A_boost)
gds_ID_cas_eff = gds_ID_cas ./ (1 + A_boost);

% 4) compute Av0 for each combination (vectorized)
% we need gds_load: p-load gds/ID at gm_ID_load for each L_cas (use pch)
gds_ID_load = zeros(size(gds_ID_cas));
for k=1:length(L_cas)
    gds_ID_load(:,k) = lookup(pch,'GDS_ID','GM_ID',gm_ID_load,'L',L_cas(k)); % p-type for load
end

% compute Av (for each L1 and L_cas)
Av = zeros(length(L1), length(L_cas));
for k=1:length(L_cas)
    Av(:,k) = gm_ID1_vec ./ ( gds_ID1_vec + gds_ID_cas_eff(:,k) + gds_ID_load(:,k) );
end

% 5) find best L1 for max |Av| for each L_cas
[Av_max, idx_L1] = max(Av,[],1);    % Av_max(k) corresponds to L1(idx_L1(k))

% 6) de-normalize (compute ID, W1, Wcas, Wload) and include self-loading iteratively
fprintf('\n--- Cascode Amplifier Design Summary ---\n');
fprintf('Target fu: %.1f GHz, Load CL: %.1f pF, Booster Gain: %.0f\n', fu/1e9, CL*1e12, A_boost);
fprintf('-------------------------------------------------\n');

Cself = 0;
for k = 1:length(L_cas)
    L1_best = L1(idx_L1(k));
    gm_ID1 = gm_ID1_vec(idx_L1(k));
    gds_ID1 = gds_ID1_vec(idx_L1(k));
    gds_IDcas = gds_ID_cas(:,k);
    gds_IDload = gds_ID_load(:,k);
    
    % iterative loop for Cself (same idea as Example 3.7/3.8)
    Cself = 0;
    for iter = 1:5
        gm_req = 2*pi*fu*(CL + Cself);   % required gm1 accounting for self-loading.
        ID = gm_req / gm_ID1;            % ID of M1
        JD1 = diag( lookup(nch,'ID_W','GM_ID',gm_ID1,'L',L1_best) );   % current density
        W1 = ID ./ JD1;
        % cascode de-normalization (assume same ID through cascode)
        JD_cas = lookup(nch,'ID_W','GM_ID',gm_ID_cas,'L',L_cas(k));
        W_cas = ID ./ JD_cas;
        % p-load de-normalization (ID flows in load too)
        JD_load = lookup(pch,'ID_W','GM_ID',gm_ID_load,'L',L_cas(k));
        W_load = ID ./ JD_load;
        % compute Cdd contributions per unit width (lookup) and total Cself
        Cdd1_w = diag( lookup(nch,'CDD_W','GM_ID',gm_ID1,'L',L1_best) );
        Cdd_cas_w = lookup(nch,'CDD_W','GM_ID',gm_ID_cas,'L',L_cas(k));
        Cdd_load_w = lookup(pch,'CDD_W','GM_ID',gm_ID_load,'L',L_cas(k));
        Cdd1 = W1 .* Cdd1_w;
        Cdd_cas = W_cas .* Cdd_cas_w;
        Cdd_load = W_load .* Cdd_load_w;
        Cself = Cdd1 + Cdd_cas + Cdd_load;
    end

    % --- Print Statements Added Here ---
    fprintf('➡️ Results for **L_cas = %.2f \u03BCm**:\n', L_cas(k)*1e6);
    fprintf('  * Max Av0: **%.2f** (%.1f dB)\n', Av_max(k), 20*log10(Av_max(k)));
    fprintf('  * Optimal L1: %.3f \u03BCm\n', L1_best*1e6);
    fprintf('  * Bias Current (ID): **%.3f \u03BCA**\n', ID*1e6);
    fprintf('  * Total C_self: **%.3f fF**\n', Cself*1e15);
    fprintf('  * Transistor Widths:\n');
    fprintf('    - W1 (Driver, N): %.1f \u03BCm\n', W1*1e6);
    fprintf('    - W_cas (Cascode, N): %.1f \u03BCm\n', W_cas*1e6);
    fprintf('    - W_load (P-Load, P): %.1f \u03BCm\n', W_load*1e6);
    fprintf('-------------------------------------------------\n');
    % -----------------------------------

    % record results per L_cas(k): Av_max(k), L1_best, W1, W_cas, W_load, ID, Cself
end
