% --- START: Av0_vs_L_and_gmID_vs_L_CASCODE.m ---
% Cascode op-amp gain evaluation using Paul Jespers' gm/ID method
% Av0 = gm1*gm3^2 / ( gm5*gds3*gds1 + gm3*gds5*gds7 )
clear; close all; clc;

% Load lookup-table structures (must contain nch and pch)
load('nch_18.mat','nch');
load('pch_18.mat','pch');

% === Design parameters ===
Lvec = 0.18:0.01:1.98;      % µm sweep
fu   = 1e9;                 % unity-gain freq target (Hz)
FO   = 10;                  % fan-out (as in Jespers Example 3.7)
gm_ID5 = 10;                % gm/ID (S/A) for pMOS cascode & current source
VSB = 0;

target_GM_CGG = 2*pi*fu*FO; % for nMOS gm/ID lookup

% === Preallocation ===
N = numel(Lvec);
gmID1 = nan(N,1); gmID3 = nan(N,1); gmID5 = nan(N,1);
gdsID1 = nan(N,1); gdsID3 = nan(N,1);
gdsID5 = nan(N,1); gdsID7 = nan(N,1);
Av0_vec = nan(N,1);

% === Sweep gate length ===
for k = 1:N
    Lk = Lvec(k);

    % --- Input device (M1) ---
    gmID1(k)  = lookup(nch,'GM_ID','GM_CGG',target_GM_CGG,'L',Lk,'VSB',VSB);
    gdsID1(k) = lookup(nch,'GDS_ID','GM_ID',gmID1(k),'L',Lk);

    % --- Cascode nMOS (M3) ---
    gmID3(k)  = gmID1(k);   % often same inversion level
    gdsID3(k) = lookup(nch,'GDS_ID','GM_ID',gmID3(k),'L',Lk);

    % --- pMOS cascode (M5) ---
    gmID5(k)  = gm_ID5;
    gdsID5(k) = lookup(pch,'GDS_ID','GM_ID',gmID5(k),'L',Lk);

    % --- pMOS current-source load (M7) ---
    gdsID7(k) = lookup(pch,'GDS_ID','GM_ID',gmID5(k),'L',Lk);

    % --- Compute normalized |Av0| (Eq. for cascode op-amp) ---
    Av0_vec(k) = (gmID1(k) .* gmID3(k).^2) ./ ...
        ( gmID5(k).*gdsID3(k).*gdsID1(k) + gmID3(k).*gdsID5(k).*gdsID7(k) );
end

% === Plot results ===
figure;
plot(Lvec, Av0_vec,'-o','LineWidth',1.6);
xlabel('L (\mum)');
ylabel('|A_{v0}| (normalized)');
title('|A_{v0}| vs L for Cascode Op-Amp (gm/ID method)');
grid on;

figure;
plot(Lvec, gmID1,'-s','LineWidth',1.6);
xlabel('L (\mum)');
ylabel('gm/ID (S/A)');
title('gm/ID (Input device M1) vs L');
grid on;

% === Save data ===
save('Av0_gmID_vs_L_CASCODE.mat','Lvec','Av0_vec',...
     'gmID1','gmID3','gmID5','gdsID1','gdsID3','gdsID5','gdsID7');
csvwrite('Av0_vs_L_CASCODE.csv',[Lvec(:) Av0_vec(:)]);
csvwrite('gmID_vs_L_CASCODE.csv',[Lvec(:) gmID1(:)]);

fprintf('✅ Cascode computation done. Results saved to MAT and CSV files.\n');
% --- END ---
