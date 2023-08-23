clc;
Ca= [0.5, 3, 1, 2, 6, 10, 8, 4];
Ca0 = [2, 5, 6, 6, 11, 14, 16, 24];
t = [30, 1, 50, 8, 4, 20, 20, 4];
vol_rate = 0.1;
Ca_initial = 10;
Xaf = 0.9;
Ca_out = Ca_initial*(1-Xaf);

Ra = t./(Ca0 - Ca); % 1/(-ra)
Xa = (Ca0 - Ca)./Ca0;

% plot Ra vs Ca
pp = spline(Ca, Ra); % represents spline function.
xx = linspace(min(Ca), max(Ca), 1000); % equally spaced x-values.
f = ppval(pp, xx); % equally spaced function values corresponding to values in xx array.

plot(Ca, Ra, '*', xx, f,'m','LineWidth',2);
legend('Data', 'Spline Fit');
grid on;
xlabel('C_A_f');
ylabel('-1/r_A');
title('Spline Interpolation of -1/rA vs CA_f');

% Calculating volumes for different reactors
V_pfr = vol_rate * integral(@(x) ppval(pp, x), Ca_out, Ca_initial);
V_cstr = vol_rate * (Ca_initial - Ca_out) * ppval(pp, Ca_out);

% MFR and PFR
Ca_min = fminbnd(@(x) ppval(pp, x), Ca_out, Ca_initial);
rate_min = ppval(pp, Ca_min);
V_pfr_mfr = vol_rate * (integral(@(x) ppval(pp, x), Ca_out, Ca_min) + (Ca_initial - Ca_min) * rate_min);

% Two CSTR
xval = Ca_out:0.01:Ca_initial;
area = @(c) (Ca_initial - c).*ppval(pp, c) + (c- Ca_out).*ppval(pp,Ca_out);
[min_area,Ca_2] = min(area(xval));
V_2cstr = min_area* vol_rate;


% recycle PFR
for Xai = 0:0.001:Xaf
    Cai = Ca_initial * (1 - Xai);
    lhs = ppval(pp, Cai) * (Xaf - Xai) * Ca_initial;
    rhs = integral(@(x) ppval(pp, x), Ca_out, Cai);
    if abs(lhs - rhs) < 0.1
        recycle_ratio = (Ca_initial - Cai) / (Cai - Ca_out);
        V_recycle = (Ca_initial - Ca_out) * ppval(pp, Cai) * vol_rate;
        volumetric_flow_recycle = vol_rate * recycle_ratio;
        break;
    end
end

% Displaying the results 
fprintf('Volume of Single PFR is %f\n', V_pfr);
fprintf('Volume of Single CSTR is %f\n', V_cstr);
fprintf('Volume of two stirred tank is %f\n', V_2cstr);
fprintf('Volume of combination of a PFR and a MFR is %f\n', V_pfr_mfr);
fprintf('Volume of a PFR with recycle is %f\n', V_recycle);
fprintf('Optimum Recycle ratio of PFR is %f\n', recycle_ratio);
fprintf('Optimum inlet Concentration(Cai- according to notes) is %f\n', Cai);
