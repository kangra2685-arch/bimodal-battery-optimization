close all ; clear ; clc ;

% --- 1. 고정 파라미터 설정 ---
R = 8.31446 ; % C*V/(K*mol)                     쿨롱*볼트/(켈빈*몰)
T = 298.15 ; % K                                켈빈
F = 96485 ; % C/mol                             쿨롱/몰
d1 = 4 * 10^(-6) ; % m                          미터
d2 = 8 * 10^(-6) ; % m                          미터
density = 4.81756 * 10^6 ; % g/m^3              그램/미터^3
D = 1.0 * 10^(-14) ; % m^2/s                    미터^2/초
i0 = 1.5; % A/m^2                               암페어/미터^2
Cdl_hat = 0.05 ; % F/m^2                        패럿/미터^2
Rs = 5; % ohm                                   옴
Vocv = 4.2; % V                                 볼트
Vmin = 3.0; % V                                 볼트
t_pulse = 30; % s                               초
total_mass_mg = 10; % mg                        밀리그램
total_mass_kg = total_mass_mg * 10^(-6); % kg   킬로그램
specific_cap_F = 600; % F/g                     패럿/그램 (비용량)
total_mass = total_mass_mg * 10^(-3); % g       그램
% 1 F = 1 C/V

% --- 2. 주파수 및 비율 설정 ---
f_min = 0.001 ; f_max = 10000 ; num_points = 70 ;
frequency_vec = logspace( log10( f_min ), log10( f_max ), num_points ) ;
w_vec = 2 * pi * frequency_vec ;
ratios = 0:0.1:1; 
colors = jet(length(ratios)); 

% --- [추가] 노이즈 파라미터 설정 ---
noise_level = 0.005; % 임피던스 크기의 0.5% 수준 노이즈

% --- 3. 시뮬레이션 루프 ---
figure('Name', 'NCM811 Specific Power Analysis (With Noise)'); hold on;

for i = 1:length(ratios)
    w1_ratio = ratios(i);
    w2_ratio = 1 - w1_ratio;
    
    x = total_mass * w1_ratio; 
    y = total_mass * w2_ratio; 
    
    S1 = (x > 0) * (3 * x) / (density * (d1/2)); 
    S2 = (y > 0) * (3 * y) / (density * (d2/2));
    Stotal = S1 + S2;
    
    CL1 = x * specific_cap_F;
    CL2 = y * specific_cap_F;
    
    % --- Branch 1 & 2 로직 (기존과 동일) ---
    if x > 0
        Rct1 = (R*T) / (F*i0*S1); tau1 = (d1/2)^2 / D;
        Y1_freq = 1 ./ (Rct1 + (tau1/(CL1*sqrt(2))) .* (coth(sqrt(1j*w_vec*tau1)) ./ sqrt(1j*w_vec*tau1)));
    else, Y1_freq = 0; end
    
    if y > 0
        Rct2 = (R*T) / (F*i0*S2); tau2 = (d2/2)^2 / D;
        Y2_freq = 1 ./ (Rct2 + (tau2/(CL2*sqrt(2))) .* (coth(sqrt(1j*w_vec*tau2)) ./ sqrt(1j*w_vec*tau2)));
    else, Y2_freq = 0; end
    
    % --- 전체 임피던스 계산 ---
    Y_total = Y1_freq + Y2_freq + (1j * w_vec * Cdl_hat * Stotal);
    Z_clean = Rs + (1 ./ Y_total); 
    
    % --- [추가] 노이즈 생성 로직 ---
    % 실수부와 허수부에 각각 독립적인 가우시안 노이즈 추가
    % 노이즈 크기는 각 지점의 임피던스 절댓값(Magnitude)에 비례하도록 설정
    Z_mag = abs(Z_clean);
    noise_real = randn(size(Z_clean)) .* Z_mag * noise_level;
    noise_imag = randn(size(Z_clean)) .* Z_mag * noise_level;
    Z_noise = real(Z_clean) + noise_real + 1j*(imag(Z_clean) + noise_imag);
    
    % --- 데이터 저장 (파일 이름에 '노이즈' 추가) ---
    % ratios(i) 값을 포함한 파일명 생성
    folderName = 'EIS_Result';
    fileName = sprintf('NCM811_EIS_Ratio_%.1f_noise.csv', w1_ratio);
    if ~isfolder(folderName)
        mkdir(folderName);
        fprintf('폴더 "%s"를 생성했습니다.\n', folderName);
    end
    filePath = fullfile(folderName, fileName);

    data_to_save = [frequency_vec', real(Z_noise)', imag(Z_noise)'];
    res_table = array2table(data_to_save, 'VariableNames', {'Frequency_Hz', 'Z_real_Ohm', 'Z_imag_Ohm'});
    writetable(res_table, filePath);

    % --- EIS 플로팅 ---
    plot(real(Z_noise), -imag(Z_noise), 'o-', 'Color', colors(i,:), ...
         'MarkerSize', 3, 'LineWidth', 0.8, 'DisplayName', sprintf('Ratio %.1f (Noise)', w1_ratio));
end

fprintf('EIS 가상 데이터를 "%s"에 저장했습니다.\n', folderName);

% --- 4. 그래프 서식 설정 ---
xlabel('Z'' (\Omega)'); ylabel('-Z'''' (\Omega)');
grid on; axis equal;
title('NCM811 EIS with Artificial Noise');
legend('Location', 'northeastoutside');