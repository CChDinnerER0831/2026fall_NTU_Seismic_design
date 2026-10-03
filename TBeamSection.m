classdef TBeamSection
    properties
        FlangeWidth = 100;      % 翼板寬度
        FlangeThickness = 10;   % 翼板厚度
        TotalHeight = 200;      % 總高度
        WebThickness = 12;      % 腹板厚度
        elasticLimmit = 0.00125;   % 彈性極限 (降伏應變)
        AC                      % 形心座標
        E = 200                 %young's module
    end
    
    methods
        function obj = TBeamSection()
            obj.AC = obj.getCentroid();
        end 
        
        function AC = getCentroid(obj)
            % 使用 integral2 進行數值面積分，完全依賴 checkPoint()
            
            % 1. 設定積分邊界 (涵蓋整個截面的 Bounding Box)
            xmin = -obj.FlangeWidth;
            xmax = obj.FlangeWidth;
            ymin = 0;
            ymax = obj.TotalHeight;
            
            % 2. 定義積分函數 (Integrands)
            funArea = @(x,y) obj.checkPoint(x,y);
            funQx   = @(x,y) y .* obj.checkPoint(x,y);
            funQy   = @(x,y) x .* obj.checkPoint(x,y);
            
            % 3. 執行內建數值積分
            Area = integral2(funArea, xmin, xmax, ymin, ymax);
            Qx = integral2(funQx, xmin, xmax, ymin, ymax);
            Qy = integral2(funQy, xmin, xmax, ymin, ymax);
            
            % 4. 計算形心
            if Area > 0
                xc = Qy / Area;
                yc = Qx / Area;
            else
                xc = 0;
                yc = 0;
            end
            
            AC = [xc, yc];
        end 
        
        function e = epsilon(obj, D, x, y)
            % 計算應變分佈，加入 .* 支援矩陣運算
            phi_x = D(1);
            phi_y = D(2);
            es = D(3);
            
            ey = phi_x .* (y - obj.AC(2));
            ex = phi_y .* (x - obj.AC(1));
            e = ex + ey + es;
        end 
        
        function yield_flag = is_yield(obj, D, x, y)
            % 判斷是否降伏。加入參數 D，並使用 abs() 考慮拉壓皆會降伏
            % 使用邏輯判斷轉 double 以支援矩陣運算
            strain = obj.epsilon(D, x, y);
            yield_flag = double(abs(strain) >= obj.elasticLimmit);
        end
        
        function is_exist = checkPoint(obj, x, y)
            % checkPoint 檢查點 (x,y) 是否存在於 T 型截面內
            % 【重要修正】：將 && 改為 &，|| 改為 |，以支援 integral2 的矩陣輸入
            
            % 1. 判斷是否在翼板 (Flange) 內
            inFlange = (y >= (obj.TotalHeight - obj.FlangeThickness)) & ...
                       (y <= obj.TotalHeight) & ...
                       (x >= -obj.FlangeWidth/2) & ...
                       (x <= obj.FlangeWidth/2);
                       
            % 2. 判斷是否在腹板 (Web) 內
            inWeb = (y >= 0) & ...
                    (y < (obj.TotalHeight - obj.FlangeThickness)) & ...
                    (x >= -obj.WebThickness/2) & ...
                    (x <= obj.WebThickness/2);
                    
            % 3. 若在翼板或腹板內，則回傳 1，否則回傳 0
            is_exist = double(inFlange | inWeb);
        end

function Q = Q_matrix(obj, D)
            % Q_matrix 計算截面切線剛度矩陣 (Tangent Stiffness Matrix)
            % D 為變形參數陣列 [phi_x, phi_y, es]
            
            % 1. 建立網格座標
            dx = 0.5;
            dy = 0.5;
            x_range = -obj.FlangeWidth : dx : obj.FlangeWidth;
            y_range = 0 : dy : obj.TotalHeight;
            [X, Y] = meshgrid(x_range, y_range);
            
            % 2. 判斷各網格點的狀態
            % 存在矩陣 (1: 存在, 0: 不存在)
            exist_matrix = obj.checkPoint(X, Y); 
            
            % 降伏矩陣 (1: 降伏, 0: 未降伏)
            yield_matrix = obj.is_yield(D, X, Y); 
            
            % 彈性有效矩陣：必須「存在」且「未降伏」
            % 根據公式，當降伏時 E = 0，因此只計算未降伏區域
            active_elastic = exist_matrix .* (1 - yield_matrix); 
            
            % 3. 定義相對於形心的局部座標 (Lever arms)
            % 應變公式是基於相對形心的距離，因此積分的 x, y 也是相對於形心
            x_bar = X - obj.AC(1);
            y_bar = Y - obj.AC(2);
            
            % 4. 準備積分參數
            dA = dx * dy;
            E_eff = obj.E .* active_elastic; % 每個點的有效 E 值 (降伏處或非截面處為 0)
            
            % 5. 根據圖片公式計算 Qij 各項元素
            Q11 = sum( E_eff(:) .* y_bar(:).^2 ) * dA;
            Q22 = sum( E_eff(:) .* x_bar(:).^2 ) * dA;
            Q33 = sum( E_eff(:) ) * dA;
            
            Q12 = sum( E_eff(:) .* x_bar(:) .* y_bar(:) ) * dA;
            Q21 = Q12;
            
            Q13 = -sum( E_eff(:) .* y_bar(:) ) * dA;
            Q31 = Q13;
            
            Q23 = -sum( E_eff(:) .* x_bar(:) ) * dA;
            Q32 = Q23;
            
            % 6. 組裝並回傳 Q 矩陣
            Q = [Q11, Q12, Q13;
                 Q21, Q22, Q23;
                 Q31, Q32, Q33];
end

        function force = get_Force(obj, D)
            % get_Force 計算給定變形參數 D [phi_x, phi_y, es] 下的對應內力
            % 回傳 force = [F1(Mx); F2(My); F3(P)]
            
            % 1. 建立網格座標
            dx = 0.05;
            dy = 0.05;
            x_range = -obj.FlangeWidth : dx : obj.FlangeWidth;
            y_range = 0 : dy : obj.TotalHeight;
            [X, Y] = meshgrid(x_range, y_range);
            
            % 2. 取得網格點狀態與應變
            exist_matrix = obj.checkPoint(X, Y);     % 存在矩陣 (1:存在, 0:不存在)
            yield_matrix = obj.is_yield(D, X, Y);    % 降伏矩陣 (1:降伏, 0:未降伏)
            strain = obj.epsilon(D, X, Y);           % 截面各點應變
            
            % 3. 計算各網格點的應力 (Stress)
            % 未降伏區域應力 (Hooke's Law: stress = E * strain)
            elastic_stress = obj.E .* strain .* (1 - yield_matrix); 
            
            % 降伏區域應力 (達到降伏上限，並保持正負號)
            fy = obj.E * obj.elasticLimmit;          % 計算降伏應力 Fy
            plastic_stress = sign(strain) .* fy .* yield_matrix; 
            
            % 總真實應力 (僅計算存在於截面內的點)
            stress = (elastic_stress + plastic_stress) .* exist_matrix; 
            
            % 4. 定義相對於形心的局部力臂
            x_bar = X - obj.AC(1);
            y_bar = Y - obj.AC(2);
            
            % 5. 數值積分計算內力
            dA = dx * dy;
            
            % F1 (對應 phi_x 的彎矩)：積分 stress * y_bar
            F1 = sum( stress(:) .* y_bar(:) ) * dA;
            
            % F2 (對應 phi_y 的彎矩)：積分 stress * x_bar
            F2 = sum( stress(:) .* x_bar(:) ) * dA;
            
            % F3 (對應 es 的軸力)：積分 stress
            F3 = sum( stress(:) ) * dA;
            
            % 6. 組裝回傳
            force = [F1; F2; F3];
        end
    end
end