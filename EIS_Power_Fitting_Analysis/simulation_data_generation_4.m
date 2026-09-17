close all ; clear ; clc ;

%% 1. 설정 및 물리 상수 (참값 설정)
R = 8.31446 ; T = 298.15 ; F = 96485 ;
d1 = 4 * 10^(-6) ; d2 = 8 * 10^(-6) ; 
density = 4.81756 * 10^6 ; 
D = 1.0 * 10^(-14) ; 

i0 = 1.5; 
Cdl_hat = 0.05 ; 
Rs_true = 5; 

Vocv = 4.2; Vmin = 3.0; 
t_short = 30;   % 30초
t_long = 1800;  % 30분
total_mass_mg = 10; 
total_mass_kg = total_mass_mg * 10^(-6); 
total_mass = total_mass_mg * 10^(-3); 
specific_cap_F = 600; 

% 고정 시상수 (참값)
tau1_true = (d1/2)^2 / D; % 400s
tau2_true = (d2/2)^2 / D; % 1600s

ratios = 0.0:0.1:1.0; 
num_ratios = length(ratios);

% 결과 저장용 행렬 (참값 보관용)
true_params_all = zeros(num_ratios, 11); % [Rs, Qdl, n, Rct1, Rw1, tau1, Rct2, Rw2, tau2, P30s, P30m]

%% 2. 시뮬레이션 및 데이터 생성 루프
for i = 1:num_ratios
    w1_ratio = ratios(i);
    w2_ratio = 1 - w1_ratio;
    
    x = total_mass * w1_ratio; 
    y = total_mass * w2_ratio; 
    
    S1 = (x > 0) * (3 * x) / (density * (d1/2)); 
    S2 = (y > 0) * (3 * y) / (density * (d2/2));
    Stotal = S1 + S2;
    
    CL1 = x * specific_cap_F;
    CL2 = y * specific_cap_F;
    
    % --- Branch 1 (Small) ---
    if x > 0
        Rct1 = (R*T) / (F*i0*S1);
        Rw1 = tau1_true / (CL1 * sqrt(2)); 
        % 30s 및 30min 저항 기여분 (가상데이터 생성 로직 기반)
        Rd1_30s = Rw1 * sqrt((4 * t_short) / (pi * tau1_true));
        Rd1_30m = Rw1 * (t_long/tau1_true + 1/3);
    else
        Rct1 = 1e10; Rw1 = 0; Rd1_30s = 1e10; Rd1_30m = 1e10; % Inf 대응
    end
    
    % --- Branch 2 (Large) ---
    if y > 0
        Rct2 = (R*T) / (F*i0*S2);
        Rw2 = tau2_true / (CL2 * sqrt(2));
        Rd2_30s = Rw2 * sqrt((4 * t_short) / (pi * tau2_true));
        Rd2_30m = Rw2 * (t_long/tau2_true + 1/3);
    else
        Rct2 = 1e10; Rw2 = 0; Rd2_30s = 1e10; Rd2_30m = 1e10;
    end
    
    % --- 전체 병렬 저항 계산 (합성 저항) ---
    % 30초 저항 및 출력
    R_path1_30s = Rct1 + Rd1_30s;
    R_path2_30s = Rct2 + Rd2_30s;
    R_tot_30s = Rs_true + 1 / (1/R_path1_30s + 1/R_path2_30s);
    P30s_true = (Vmin * (Vocv - Vmin)) / (R_tot_30s * 1e-5); % 기존 단위 보정 반영
    
    % 30분 저항 및 출력
    R_path1_30m = Rct1 + Rd1_30m;
    R_path2_30m = Rct2 + Rd2_30m;
    R_tot_30m = Rs_true + 1 / (1/R_path1_30m + 1/R_path2_30m);
    P30m_true = (Vmin * (Vocv - Vmin)) / (R_tot_30m * 1e-5);

    % --- 참값 저장 ---
    % [Rs, Qdl, n, Rct1, Rw1, tau1, Rct2, Rw2, tau2, P30s, P30m]
    % (가상데이터에서 n은 1.0, Qdl은 Cdl_hat*Stotal로 설정됨)
    true_params_all(i, :) = [Rs_true, Cdl_hat*Stotal, 1.0, Rct1, Rw1, tau1_true, Rct2, Rw2, tau2_true, P30s_true, P30m_true];
    
    % --- 가상 데이터 파일 저장 (피팅 코드에서 읽을 파일) ---
    % (실제 저장 로직은 생략하거나 아래와 같이 테이블로 생성 가능)
    % csv_data = table(frequency_vec', real(Z_final)', imag(Z_final)', ...);
    % writetable(csv_data, sprintf('EIS_Data_Ratio_%.1f.csv', w1_ratio));
end

%% 3. 결과 리포트 출력 (가상 데이터 기반 '참값' 요약)
fprintf('\n==========================================================================================\n');
fprintf('   [가상 데이터 생성에 사용된 파라미터 참값(Original True Values) 리포트]   \n');
fprintf('==========================================================================================\n');

param_names = {'Rs', 'Qdl', 'n', 'Rct1', 'Rw1', 'tau1', 'Rct2', 'Rw2', 'tau2', 'P_30s', 'P_30min'};
True_Table = array2table(true_params_all, ...
    'VariableNames', param_names, ...
    'RowNames', string(ratios));

% Inf 또는 너무 큰 값 가독성 처리
disp(True_Table);

%% 4. 비교를 위한 간단한 시각화
figure('Name', 'True Parameter vs Ratio');
subplot(1,2,1);
plot(ratios, true_params_all(:,4), 'o-', 'DisplayName', 'Rct1 (Small)'); hold on;
plot(ratios, true_params_all(:,7), 's-', 'DisplayName', 'Rct2 (Large)');
set(gca, 'YScale', 'log'); grid on;
title('True Charge Transfer Resistance'); xlabel('Ratio'); ylabel('\Omega');
legend;

subplot(1,2,2);
plot(ratios, true_params_all(:,10), 'b-o', 'LineWidth', 2, 'DisplayName', '30s Power'); hold on;
plot(ratios, true_params_all(:,11), 'r-s', 'LineWidth', 2, 'DisplayName', '30m Power');
grid on; title('True Specific Power Trend'); xlabel('Ratio'); ylabel('W/kg');
legend;

fprintf('\n* 주의: Ratio 0.0 또는 1.0에서 나타나는 매우 큰 값(1e10)은 해당 경로가 존재하지 않음을 의미합니다.\n');