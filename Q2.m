clear;close all;clc
% 1. 建立 T 型梁物件並設定正確的物理參數 (單位: N, mm, MPa)
T = TBeamSection();
T.E = 200000;               % 彈性模數 200 GPa = 200,000 MPa (N/mm^2)
T.elasticLimmit = 0.00125;  % 降伏應變 = 250 MPa / 200,000 MPa

% 2. 設定目標軸力 (壓力為負)
P_target = -164000;         % -164 kN 轉換為 N

% 3. 設定位移控制的曲率 (phi_x) 增量步
% 曲率從 0 逐步增加至 0.00006 rad/mm，切分 60 步
phi_steps = linspace(0, 0.00006, 60);
M_values = zeros(1, length(phi_steps)); % 儲存對應的彎矩

% 初始化變形向量 [phi_x, phi_y, es]
D = [0, 0, 0];

% 4. 開始位移控制迴圈
for i = 1:length(phi_steps)
    % 固定當前步的曲率
    D(1) = phi_steps(i); 
    D(2) = 0; % 假設無雙軸彎曲
    
    iter = 0;
    max_iter = 100;
    tolerance = 10; % 容許力誤差 10 N
    
    % 針對軸力進行 1D Newton-Raphson 疊代 (尋找正確的中性軸位置 es)
    while iter < max_iter
        % 計算當前應變狀態下的內部力
        Fi = T.get_Force(D);
        P_current = Fi(3); 
        
        % 計算軸力殘餘值 (Residual)
        residual = P_target - P_current;
        
        % 若達到容許誤差內，跳出疊代
        if abs(residual) < tolerance
            break;
        end
        
        % 計算剛度矩陣，取出 dP/des (即 Q33)
        Q = T.Q_matrix(D);
        K33 = Q(3, 3);
        
        % 防呆：若全截面降伏，K33 會趨近 0
        if K33 < 1e-5
            warning('截面於步數 %d 失去軸向剛度 (已完全塑性)', i);
            break;
        end
        
        % 更新軸向應變 es (D(3))
        D(3) = D(3) + (residual / K33);
        
        iter = iter + 1;
    end
    
    % 將收斂後的內部彎矩 Mx 儲存起來 (N-mm 轉換為 kN-m 方便繪圖)
    M_values(i) = Fi(1) / 1e6;
end

% 5. 繪製 M-phi 曲線
figure;
plot(phi_steps, M_values, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
grid on;
xlabel('曲率 \phi_x (rad/mm)', 'FontSize', 12);
ylabel('彎矩 M_x (kN-m)', 'FontSize', 12);
title('T型梁 M-\phi 關係曲線 (軸力 P = -164 kN)', 'FontSize', 14);

% 標示出最大彎矩 (塑性彎矩容量)
max_M = max(M_values);
yline(max_M, '--r', sprintf('M_{max} \\approx %.2f kN-m', max_M), 'LabelHorizontalAlignment', 'left');