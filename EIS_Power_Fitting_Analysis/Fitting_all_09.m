clear; clc; close all;

%% 1. 설정 및 물리 상수 계산
ratios = 0.0:0.1:1.0; 
file_pattern = 'EIS_Data_Ratio_%.1f.csv';
t_short = 30;    
t_long = 1800;   
V_min = 2.5;  
V_ocv = 3.7;  

% [입력받은 물리 상수 적용]
D = 1.0e-14; % 확산계수 (m^2/s)
R1 = 2.0e-6; % 소립자 반경 (4um 지름 -> 2um 반경)
R2 = 4.0e-6; % 대립자 반경 (8um 지름 -> 4um 반경)

tau1_const = (R1^2) / D; % 400s
tau2_const = (R2^2) / D; % 1600s

fprintf('계산된 물리 상수: tau1(소립)=%.1fs, tau2(대립)=%.1fs\n', tau1_const, tau2_const);

results_R2 = zeros(length(ratios), 1);
specific_power_30s = zeros(length(ratios), 1);
specific_power_30min = zeros(length(ratios), 1);
fit_params_all = zeros(length(ratios), 9); 

figure('Name', 'EIS Fitting with Physical Tau (400s/1600s)', 'NumberTitle', 'off', 'Position', [50, 50, 1400, 700]);

%% 2. 루프: 데이터 로드 및 피팅
for i = 1:length(ratios)
    ratio = ratios(i);
    filename = sprintf(file_pattern, ratio);
    
    if ~exist(filename, 'file'), continue; end
    data = readtable(filename);
    freq = data.Frequency_Hz;
    z_real = data.Real_Z;
    z_imag = data.Imag_Z;
    z_exp = z_real + 1j*z_imag;
    omega = 2 * pi * freq;
    
    % [모델] p = [Rs, Qdl, n, Rct1, Rw1, Rct2, Rw2] (tau는 상수로 고정)
    fit_func = @(p, w) p(1) + 1 ./ (p(2)*(1j*w).^p(3) + ...
        (1./(p(4) + p(5)*coth(sqrt(1j*w*tau1_const))./sqrt(1j*w*tau1_const)) + ...
         1./(p(6) + p(7)*coth(sqrt(1j*w*tau2_const))./sqrt(1j*w*tau2_const))));

    weight = 1 ./ abs(z_exp); 
    obj_func = @(p, w) [real(fit_func(p, w)).*weight, imag(fit_func(p, w)).*weight];
    y_weighted = [z_real.*weight, z_imag.*weight];

    % 초기값 설정 (물리적 일관성을 위해 비율에 따라 Rct1 하향 유도)
    R_est = max(z_real) - min(z_real);
    p0 = [min(z_real), 1e-4, 0.85, R_est*(1-ratio*0.5), R_est*0.5, R_est*(0.5+ratio), R_est*0.5]; 
    lb = [0.1, 1e-8, 0.5, 0.001, 0.001, 0.001, 0.001];
    ub = [100, 0.5, 1.0, 3000, 3000, 3000, 3000];

    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'FunctionTolerance', 1e-14, 'StepTolerance', 1e-14, 'MaxFunctionEvaluations', 10000);
    
    [p_res, ~] = lsqcurvefit(obj_func, p0, omega, y_weighted, lb, ub, opts);
    
    % 파라미터 저장 (출력용 9개 구성)
    p_full = [p_res(1:5), tau1_const, p_res(6:7), tau2_const];
    fit_params_all(i, :) = p_full;
    
    z_fit = fit_func(p_res, omega);
    results_R2(i) = 1 - (sum(abs(z_exp - z_fit).^2) / sum(abs(z_exp - mean(z_exp)).^2));

    %% 3. 비출력 계산
    % 30초 출력 (Semi-infinite)
    Rd1_30s = p_full(5) * sqrt((4 * t_short) / (pi * tau1_const));
    Rd2_30s = p_full(8) * sqrt((4 * t_short) / (pi * tau2_const));
    R_tot_30s = p_full(1) + ((p_full(4)+Rd1_30s)*(p_full(7)+Rd2_30s)) / ((p_full(4)+Rd1_30s)+(p_full(7)+Rd2_30s));
    specific_power_30s(i) = (V_min * (V_ocv - V_min)) / (R_tot_30s * 1e-5);

    % 30분 출력 (Long-term diffusion)
    Rd1_30m = p_full(5) * (t_long/tau1_const + 1/3);
    Rd2_30m = p_full(8) * (t_long/tau2_const + 1/3);
    R_tot_30min = p_full(1) + ((p_full(4)+Rd1_30m)*(p_full(7)+Rd2_30m)) / ((p_full(4)+Rd1_30m)+(p_full(7)+Rd2_30m));
    specific_power_30min(i) = (V_min * (V_ocv - V_min)) / (R_tot_30min * 1e-5);

    %% 4. 시각화
    subplot(3, 4, i);
    plot(z_real, -z_imag, 'bo', 'MarkerSize', 2); hold on;
    plot(real(z_fit), -imag(z_fit), 'r-', 'LineWidth', 1.2);
    title(sprintf('Ratio %.1f (R^2: %.4f)', ratio, results_R2(i)));
    grid on; axis tight;
end

%% 5. 통합 트렌드 분석
figure('Name', 'Physical Analysis Result', 'Position', [100, 100, 1100, 500]);
subplot(1, 2, 1);
yyaxis left; plot(ratios, specific_power_30s, 'o-b', 'LineWidth', 2); ylabel('30s Power');
yyaxis right; plot(ratios, specific_power_30min, 's--r', 'LineWidth', 2); ylabel('30min Power');
grid on; title('Specific Power vs. Ratio'); xlabel('Small Particle Ratio');
legend('30s Short-term', '30min Long-term');

subplot(1, 2, 2);
plot(ratios, fit_params_all(:, 4), 'd-k', 'LineWidth', 2, 'MarkerFaceColor', 'g'); hold on;
plot(ratios, fit_params_all(:, 7), 'v-m', 'LineWidth', 2, 'MarkerFaceColor', 'm');
grid on; title('Rct Trends (Fixed Tau)'); xlabel('Ratio');
ylabel('Resistance (\Omega)'); legend('Rct1 (Small)', 'Rct2 (Large)');

%% 6. 결과 리포트
res_table = array2table([fit_params_all, results_R2, specific_power_30s, specific_power_30min], ...
    'VariableNames', {'Rs', 'Qdl', 'n', 'Rct1', 'Rw1', 'tau1', 'Rct2', 'Rw2', 'tau2', 'R2', 'P_30s', 'P_30min'}, ...
    'RowNames', string(ratios));
disp(res_table);