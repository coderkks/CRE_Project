% Define data
ca0 = [2 5 6 6 11 14 16 24];
ca_in = 10;
ca_out = 1;
caf = [0.5 3 1 2 6 10 8 4];
rate = [20 0.5 10 2 0.8 5 2.5 0.2];
t = [30 1 50 8 4 20 20 4];

% Calculate parameters
xa = (ca0 - caf) ./ ca0;
r_a = (ca0 - caf) ./ t;
rate_inv = 1 ./ r_a;

% Fit spline curve to data
pp = spline(caf, rate_inv);

% Evaluate spline curve at a fine grid of points
xgrid = linspace(min(caf), max(caf), 1000);
f = ppval(pp, xgrid);

% Plot data points and spline curve
figure;
scatter(caf, rate_inv, 'filled');
hold on;
plot(xgrid, f, 'r-', 'LineWidth', 2);
grid on;
xlabel('CAf');
ylabel('-1/rA');
title('Spline Interpolation of -1/rA vs CA');
legend('Data', 'Spline Curve');

% Calculate volumes for different reactors
V0 = 0.1;
V_pfr = V0 * integral(@(x) ppval(pp, x), ca_out, ca_in);
V_cstr = V0 * (ca_in - ca_out) * ppval(pp, ca_out);
ca_min = fminbnd(@(x) ppval(pp, x), ca_out, ca_in);
rate_min = ppval(pp, ca_min);
V_pfr_mfr = V0 * (integral(@(x) ppval(pp, x), ca_out, ca_min) + (ca_in - ca_min) * rate_min);
xa_f = (ca_in - ca_out) / ca_in;
for xai = 0:0.001:xa_f
    ca_opt = ca_in * (1 - xai);
    lhs = ppval(pp, ca_in * (1 - xai)) * (xa_f - xai) * ca_in;
    rhs = integral(@(x) ppval(pp, x), ca_in * (1 - xa_f), ca_opt);
    if abs(lhs - rhs) < 0.1
        recycle_ratio = (ca_in - ca_opt) / (ca_opt - ca_out);
        V_recycle = (ca_in - ca_out) * ppval(pp, ca_opt) * V0;
        volumetric_flow_recycle = V0 * recycle_ratio;
        break;
    end
end

% Display results
fprintf('Volume of Single PFR is %f\n', V_pfr);
fprintf('Volume of Single CSTR is %f\n', V_cstr);
fprintf('Volume of two stirred tank is %f\n', V_pfr_mfr);
fprintf('Volume of combination of a PFR and a MFR is %f\n', V_pfr_mfr);
fprintf('Volume of a PFR with recycle is %f\n', V_recycle);
fprintf('Optimum Recycle ratio of PFR is %f\n', recycle_ratio);
fprintf('Optimum inlet Concentration is %f\n', ca_opt);