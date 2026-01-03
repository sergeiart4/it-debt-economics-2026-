%% СЦЕНАРИЙ 3: «ТЕХНОЛОГИЧЕСКИЙ РЫВОК 2026»
clear all; close all; clc;

%% ПАРАМЕТРЫ
R = 10; k = 0.3; gamma0 = 0.02; gamma_new = 0.005;
T = 120; dt = 0.05; t = 0:dt:T; N = length(t);
t_rush = 60; rush_duration = 10; rush_cost = 0.3; rush_D_reduction = 0.7;

%% РАСЧЕТ
V = zeros(N,1); D = zeros(N,1); gamma = gamma0 * ones(N,1);
V(1) = 0.1; D(1) = 0;

for i = 2:N
    if abs(t(i) - t_rush) < dt/2
        D(i) = D(i-1) * (1 - rush_D_reduction);
        gamma(i) = gamma_new;
        V(i) = V(i-1) * (1 - rush_cost);
        continue;
    end
    
    is_rush = (t(i) >= t_rush && t(i) < t_rush + rush_duration);
    
    if ~is_rush
        gamma(i) = gamma(i-1);
        dDdt = k * R - gamma(i) * D(i-1)^1.2;
        
        if t(i) > t_rush + rush_duration
            eff = 1.3 - 0.3 * (1 - 2./(1 + exp(2*0.03*D(i-1))));
        else
            eff = 1 - 0.3 * (1 - 2./(1 + exp(2*0.1*D(i-1))));
        end
        dVdt = R * eff;
    else
        dDdt = 0.05 * k * R;
        dVdt = 0.3 * R;
    end
    
    D(i) = D(i-1) + dt * dDdt;
    V(i) = V(i-1) + dt * dVdt;
end

%% СГЛАЖИВАНИЕ ДАННЫХ (заменяем smooth)
window = 10;
eff_smooth = zeros(N,1);
for i = 1:N
    start_idx = max(1, i-floor(window/2));
    end_idx = min(N, i+floor(window/2));
    eff_smooth(i) = mean(V(start_idx:end_idx));
end
eff_smooth = eff_smooth / max(eff_smooth) * 0.8 + 0.2;

ratio = D ./ (V + 0.1);

%% ВИЗУАЛИЗАЦИЯ
figure('Position', [100, 100, 1600, 900], 'Color', [0.05 0.05 0.08]);

% Основной 3D график
ax_main = subplot(2,3,[1 2 4 5]);
set(ax_main, 'Color', [0.08 0.08 0.12], 'GridColor', [0.3 0.3 0.4]);
hold on;

% Фоновая сетка
[Xg, Yg] = meshgrid(linspace(0,T,15), linspace(0,max(V),10));
Zg = zeros(size(Xg));
surf(Xg, Yg, Zg, 'FaceColor', [0.2 0.2 0.3], 'EdgeColor', [0.4 0.4 0.5], 'FaceAlpha', 0.1);

% 3D траектории
h_traj = plot3(t, V, D, 'w-', 'LineWidth', 4);
V_no_rush = R * t .* (1 - 0.3 * (1 - 2./(1 + exp(2*0.1*R*t*0.1))));
h_compare = plot3(t, V_no_rush, zeros(size(t)), 'g--', 'LineWidth', 2);

% Зона рывка
[x_rush, y_rush] = meshgrid([t_rush-2, t_rush+rush_duration+2], [0, max(V)*1.1]);
z_rush = zeros(2);
surf(x_rush, y_rush, z_rush, 'FaceColor', [1 0.8 0], 'EdgeColor', 'none', 'FaceAlpha', 0.08);

% Линия рывка
plot3([t_rush t_rush], [0 max(V)*1.1], [0 0], 'y--', 'LineWidth', 3);

% Точка
h_point = plot3(t(1), V(1), D(1), 'o', 'MarkerSize', 25, ...
               'MarkerFaceColor', [1 0.3 0.3], 'MarkerEdgeColor', [1 1 1], 'LineWidth', 2);

% Векторы на поверхности
[xv, yv] = meshgrid(10:15:T-10, 20:30:max(V)-10);
for xi = 1:size(xv,1)
    for yi = 1:size(xv,2)
        tx = xv(xi,yi); vy = yv(xi,yi);
        if tx < t_rush
            dz = 0.02 * tx;
            color = [0.8 0.3 0.3];
        else
            dz = 0.005 * (tx - t_rush);
            color = [0.3 0.8 0.3];
        end
        plot3([tx tx], [vy vy], [0 dz], 'Color', color, 'LineWidth', 1.5);
    end
end

% Настройка
view(50, 30); grid on; axis tight;
xlabel('Время', 'FontSize', 13, 'Color', 'w', 'FontWeight', 'bold');
ylabel('Ценность V', 'FontSize', 13, 'Color', 'w', 'FontWeight', 'bold');
zlabel('Долг D', 'FontSize', 13, 'Color', 'w', 'FontWeight', 'bold');
title('«ТЕХНОЛОГИЧЕСКИЙ РЫВОК 2026»', 'FontSize', 16, 'Color', [1 1 0.9], 'FontWeight', 'bold');

light('Position', [T, max(V), max(D)], 'Style', 'infinite');
lighting gouraud;

%% ДОПОЛНИТЕЛЬНЫЕ ГРАФИКИ
% 1. График эффективности
ax1 = subplot(2,3,3);
set(ax1, 'Color', [0.08 0.08 0.12], 'XColor', 'w', 'YColor', 'w');

h_eff = plot(t, eff_smooth, 'c-', 'LineWidth', 2.5);
hold on;
plot([t_rush t_rush], [0 1.2], 'y--', 'LineWidth', 2);
plot([t_rush+rush_duration t_rush+rush_duration], [0 1.2], 'y--', 'LineWidth', 1);

grid on; xlim([0 T]); ylim([0.4 1.2]);
xlabel('Время', 'FontSize', 10, 'Color', 'w');
ylabel('Эффективность', 'FontSize', 10, 'Color', 'w');
title('Эффективность команды', 'FontSize', 11, 'Color', 'w', 'FontWeight', 'bold');

% 2. График отношения D/V
ax2 = subplot(2,3,6);
set(ax2, 'Color', [0.08 0.08 0.12], 'XColor', 'w', 'YColor', 'w');

h_ratio = plot(t, ratio, 'm-', 'LineWidth', 2.5);
hold on;
plot([t_rush t_rush], [0 max(ratio)*1.2], 'y--', 'LineWidth', 2);
plot([0 T], [0.5 0.5], 'r--', 'LineWidth', 1.5);

grid on; xlim([0 T]); ylim([0 max(ratio)*1.2]);
xlabel('Время', 'FontSize', 10, 'Color', 'w');
ylabel('D/V', 'FontSize', 10, 'Color', 'w');
title('Отношение долга к ценности', 'FontSize', 11, 'Color', 'w', 'FontWeight', 'bold');

%% АНИМАЦИЯ
fprintf('Запуск анимации...\n');
pause(2);

% Точка окупаемости
payback_idx = [];
for i = 1:N
    if t(i) > t_rush + rush_duration && V(i) > V_no_rush(i)
        payback_idx = i;
        break;
    end
end

if ~isempty(payback_idx)
    payback_t = t(payback_idx);
    payback_V = V(payback_idx);
    payback_D = D(payback_idx);
end

skip = 4;
explosion_drawn = false;
rush_idx = find(t >= t_rush, 1);

for i = 1:skip:N
    % Обновление траектории
    set(h_traj, 'XData', t(1:i), 'YData', V(1:i), 'ZData', D(1:i));
    
    % Пульсация точки
    if t(i) < t_rush
        pulse = 25 + 8*sin(i/10);
        color = [1 0.3 0.3];
    elseif t(i) < t_rush + rush_duration
        pulse = 22 + 10*sin(i/6);
        color = [1 0.8 0];
    else
        pulse = 28 + 5*sin(i/15);
        color = [0.3 0.8 0.3];
    end
    
    set(h_point, 'XData', t(i), 'YData', V(i), 'ZData', D(i), ...
                'MarkerSize', pulse, 'MarkerFaceColor', color);
    
    % Эффект взрыва
    if t(i) >= t_rush && ~explosion_drawn
        % Большой взрыв
        n_spheres = 20;
        for s = 1:n_spheres
            [xs, ys, zs] = sphere(6);
            scale = 1 + rand()*3;
            xs = xs * scale + t_rush + (rand()-0.5)*8;
            ys = ys * scale + V(rush_idx) + (rand()-0.5)*max(V)/10;
            zs = zs * scale + D(rush_idx)*(1-rush_D_reduction)/2;
            
            col = [1, 0.5+rand()*0.5, 0];
        
                
        end
        
        % Искры
        for spark = 1:50
            ang = rand()*2*pi;
            dist = rand()*5;
            x_spark = t_rush + cos(ang)*dist;
            y_spark = V(rush_idx) + sin(ang)*dist;
            z_spark = rand()*max(D)/2;
            plot3(ax_main, [t_rush x_spark], [V(rush_idx) y_spark], ...
                  [D(rush_idx)*(1-rush_D_reduction) z_spark], ...
                  'Color', [1 0.8 0], 'LineWidth', 1, 'LineStyle', '-');
        end
        
        % Текст взрыва
        text(ax_main, t_rush, max(V)*0.7, max(D)*0.8, 'ВЗРЫВ АРХИТЕКТУРЫ!', ...
             'FontSize', 18, 'Color', [1 0.8 0], 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center', 'Rotation', 15);
        
        explosion_drawn = true;
    end
    
    % Динамическая камера
    if t(i) < t_rush
        view_az = 50 + 10*sin(i/80);
        view_el = 30 + 5*cos(i/60);
    elseif t(i) < t_rush + rush_duration
        shake = 5*sin(i/5);
        view_az = 50 + shake;
        view_el = 30 + shake*0.5;
    else
        view_az = 50 + 8*sin(i/40);
        view_el = 30 + 3*cos(i/50);
    end
    
    if ~isempty(payback_idx) && i == payback_idx
        view_az = 70;
        view_el = 20;
        campos([payback_t*1.5, payback_V*0.7, max(D)*2]);
    else
        view(view_az, view_el);
    end
    
    % Обновление графиков
    set(h_eff, 'XData', t(1:i), 'YData', eff_smooth(1:i));
    set(h_ratio, 'XData', t(1:i), 'YData', ratio(1:i));
    
    % Эффект окупаемости
    if ~isempty(payback_idx) && abs(t(i) - payback_t) < dt*2
        % Круг окупаемости
        theta = linspace(0, 2*pi, 50);
        r = 2;
        x_circle = payback_t + r*cos(theta);
        y_circle = payback_V + r*sin(theta)*max(V)/50;
        z_circle = payback_D + zeros(size(theta));
        
        plot3(ax_main, x_circle, y_circle, z_circle, ...
              'Color', [0 1 0.5], 'LineWidth', 3, 'LineStyle', '-');
        
        text(ax_main, payback_t, payback_V*1.15, payback_D, '✓ ОКУПАЕМОСТЬ', ...
             'FontSize', 14, 'Color', [0 1 0.5], 'FontWeight', 'bold');
    end
    
    % Эффект ускорения
    if t(i) > t_rush + rush_duration && mod(i,20)==0
        dV = diff(V(max(1,i-10):i));
        if ~isempty(dV) && mean(dV) > R*1.1
            scatter3(ax_main, t(i), V(i), D(i), 50, ...
                    'o', 'MarkerFaceColor', [0.3 0.8 1], ...
                    'MarkerEdgeColor', 'none');
        end
    end
    
    drawnow;
    
    % Скорость анимации
    if t(i) < t_rush
        pause(0.025);
    elseif t(i) < t_rush + rush_duration
        pause(0.015);
    else
        pause(0.01);
    end
end

%% ФИНАЛЬНЫЙ КАДР
figure('Position', [200, 200, 1200, 500], 'Color', [0.05 0.05 0.08]);

% Сравнение
subplot(1,2,1);
hold on;

if ~isempty(payback_idx)
    idx_start = find(t >= t_rush, 1);
    idx_end = payback_idx;
    
    t_fill = t(idx_start:idx_end);
    V_fill = V(idx_start:idx_end);
    Vnr_fill = zeros(size(t_fill));
    for j = 1:length(t_fill)
        [~, closest_idx] = min(abs(t - t_fill(j)));
        Vnr_fill(j) = V_no_rush(closest_idx);
    end
    
   
         
end

plot(t, V, 'b-', 'LineWidth', 4);
plot(t, V_no_rush, 'g--', 'LineWidth', 2.5);
plot([t_rush t_rush], [0 max(V)*1.1], 'y--', 'LineWidth', 3);

if ~isempty(payback_idx)
    plot([payback_t payback_t], [0 max(V)*1.1], 'm-', 'LineWidth', 2.5);
    scatter(payback_t, payback_V, 120, 'mo', 'filled', 'LineWidth', 2);
    text(payback_t, payback_V*0.9, sprintf('t=%.0f', payback_t), ...
         'Color', 'm', 'FontSize', 11, 'FontWeight', 'bold');
end

scatter(t_rush, V(rush_idx), 100, 'yo', 'filled', 'LineWidth', 2);

grid on; xlabel('Время'); ylabel('Ценность');
title('Точка окупаемости технологического рывка', 'Color', 'w', 'FontSize', 13);
set(gca, 'Color', [0.1 0.1 0.15], 'XColor', 'w', 'YColor', 'w');

% 3D финальная траектория
subplot(1,2,2);
hold on;

% Голограмма
[Xh, Yh] = meshgrid(linspace(0,T,20), linspace(0,max(V),15));
Zh = zeros(size(Xh));
for xi = 1:size(Xh,1)
    for yi = 1:size(Xh,2)
        if Xh(xi,yi) < t_rush
            Zh(xi,yi) = 0.02 * Xh(xi,yi);
        else
            Zh(xi,yi) = 0.005 * (Xh(xi,yi) - t_rush);
        end
    end
end
surf(Xh, Yh, Zh, 'FaceColor', [0.3 0.7 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.15);

plot3(t, V, D, 'w-', 'LineWidth', 4);
plot3(t, V_no_rush, zeros(size(t)), 'g--', 'LineWidth', 2);
plot3([t_rush t_rush], [0 max(V)*1.2], [0 0], 'y--', 'LineWidth', 3);

if ~isempty(payback_idx)
    scatter3(payback_t, payback_V, payback_D, 150, ...
             'mo', 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'w');
end

scatter3(t_rush, V(rush_idx), D(rush_idx), 120, ...
         'yo', 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'w');

view(60, 25); grid on; xlabel('Время'); ylabel('Ценность'); zlabel('Долг');
title('Финальная 3D траектория', 'Color', 'w', 'FontSize', 13);
set(gca, 'Color', [0.08 0.08 0.12], 'GridColor', [0.3 0.3 0.4]);

fprintf('\nТехнологический рывок завершен!\n');