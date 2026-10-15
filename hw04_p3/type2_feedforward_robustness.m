%% 前饋路徑 Type 2 系統：參數公差對斜坡追蹤誤差的影響
% -------------------------------------------------------------------------
% 系統（Fig. 4.36(b)）的追蹤誤差轉移函數：
%
%   E(s)     tau*s^2 + (1 - A*H_f)*s + A*(1 - H_r)
%   ----  =  --------------------------------------
%   R(s)              s*(tau*s + 1) + A
%
% 斜坡參考輸入 R(s) = r_o / s^2，亦即 r(t) = r_o * t。
%
% 穩態行為（暫態極點皆在左半平面，暫態會消失）：
%   e(t) → (1 - H_r)*r_o*t + r_o*(H_r/A - H_f)    當 t → ∞
%
%   * H_r = 1 且 H_f = 1/A：斜率與常數項同時為 0，e_ss = 0（理想 Type 2）
%   * 只維持 H_r = 1：斜率為 0，e_ss = r_o*(1/A - H_f)，為有限常數
%   * H_r ≠ 1：誤差以斜率 (1 - H_r)*r_o 線性成長，沒有有限穩態值
%
% 本腳本分別讓 H_r、H_f 相對理想值偏移 0%、+1%、+5%、+10%，
% 觀察 e(t)。偏移採相對理想值的比例：H = H_ideal * (1 + 偏移)。
%
% 需要 Control System Toolbox（tf、lsim）。
% -------------------------------------------------------------------------

clear; clc; close all;

%% 1. 標稱參數與理想控制器
A       = 10;       % 受控廠增益
tau     = 1;        % 受控廠時間常數
r_o     = 1;        % 斜坡斜率，r(t) = r_o * t

Hr_nom  = 1;        % 理想參考輸入增益
Hf_nom  = 1 / A;    % 理想前饋增益，使 A*H_f = 1

% 相對理想值的百分比偏移。0 代表理想參數。
deltaPct = [0, 1, 5, 10];

%% 2. 模擬時間與斜坡輸入
tEnd = 30;
dt   = 0.01;
t    = (0:dt:tEnd).';   % 欄向量，0 ~ 30 秒
r    = r_o * t;         % 斜坡輸入 r(t)

scriptDir = fileparts(mfilename('fullpath'));
cnFont    = 'Microsoft JhengHei';

%% 3. 情境一：固定 H_f 為理想值，只讓 H_r 偏移
nCase = numel(deltaPct);
tmpHr = cell(1, nCase);
for k = 1:nCase
    Hr = Hr_nom * (1 + deltaPct(k) / 100);
    Hf = Hf_nom;
    tmpHr{k} = simulateError(A, tau, Hr, Hf, r, t, r_o, deltaPct(k));
end
caseHr = [tmpHr{:}];

%% 4. 情境二：固定 H_r 為理想值，只讓 H_f 偏移
tmpHf = cell(1, nCase);
for k = 1:nCase
    Hr = Hr_nom;
    Hf = Hf_nom * (1 + deltaPct(k) / 100);
    tmpHf{k} = simulateError(A, tau, Hr, Hf, r, t, r_o, deltaPct(k));
end
caseHf = [tmpHf{:}];

%% 5. 在 Command Window 印出理論值與模擬終點值，並寫入文字檔
outTxt = fullfile(scriptDir, 'steady_state_results.txt');
fid = fopen(outTxt, 'w', 'n', 'UTF-8');
if fid < 0
    error('無法建立結果檔：%s', outTxt);
end

logf(fid, '前饋 Type 2 系統之斜坡追蹤誤差\n');
logf(fid, '標稱參數：A = %.6g，tau = %.6g，r_o = %.6g\n', A, tau, r_o);
logf(fid, '理想控制器：H_r = %.6g，H_f = 1/A = %.6g\n', Hr_nom, Hf_nom);
logf(fid, '漸近式：e(t) → (1 - H_r)*r_o*t + r_o*(H_r/A - H_f)\n');
logf(fid, '模擬終點 t = %.6g s。H_r = 1 時，e(t_end) 才是穩態誤差；\n', tEnd);
logf(fid, 'H_r ≠ 1 時誤差仍以恆定斜率成長，e(t_end) 不是穩態值。\n\n');

logf(fid, '========================================================================\n');
logf(fid, '情境一：固定 H_f = %.6g，H_r 相對理想值偏移\n', Hf_nom);
logf(fid, '========================================================================\n');
printScenario(fid, caseHr);

logf(fid, '\n========================================================================\n');
logf(fid, '情境二：固定 H_r = %.6g，H_f 相對理想值偏移\n', Hr_nom);
logf(fid, '========================================================================\n');
printScenario(fid, caseHf);

errAll = [caseHr.eSim] - [caseHr.ePred];
errAll = [errAll, [caseHf.eSim] - [caseHf.ePred]];
logf(fid, '\n模擬終點與漸近式預測的最大絕對差：%.3e\n', max(abs(errAll)));
logf(fid, '（此差值只反映尚未完全消失的暫態，不代表公式誤差。）\n');

fclose(fid);
fprintf('結果表已寫入：%s\n', outTxt);

%% 6. 繪圖：上下兩張子圖
colors = [ ...
    0.00, 0.45, 0.74;   % 藍
    0.85, 0.33, 0.10;   % 橘
    0.47, 0.67, 0.19;   % 綠
    0.49, 0.18, 0.56];  % 紫
styles = {'-', '--', '-.', ':'};

fig = figure('Name', 'Type 2 前饋強健性', 'Color', 'w', ...
    'Position', [80, 60, 960, 820]);

% ----- 圖一：H_r 偏移。曲線會分開並持續傾斜，代表誤差無界 -----
ax1 = subplot(2, 1, 1);
hold(ax1, 'on');
for k = 1:nCase
    plot(ax1, t, caseHr(k).e, 'Color', colors(k, :), ...
        'LineStyle', styles{k}, 'LineWidth', 1.6);
end
hold(ax1, 'off');
grid(ax1, 'on');
box(ax1, 'on');
xlabel(ax1, '時間 t（秒）');
ylabel(ax1, '追蹤誤差 e(t)');
title(ax1, sprintf('固定 H_f = %.4f（理想值），H_r 偏移 0%%、+1%%、+5%%、+10%%', Hf_nom));
legend(ax1, legendLabels(caseHr), 'Location', 'northeast');
styleAxes(ax1, cnFont);

% ----- 圖二：H_f 偏移。四條曲線暫態幾乎重疊，穩態只差一個小常數 -----
ax2 = subplot(2, 1, 2);
hold(ax2, 'on');
for k = 1:nCase
    plot(ax2, t, caseHf(k).e, 'Color', colors(k, :), ...
        'LineStyle', styles{k}, 'LineWidth', 1.6);
end
hold(ax2, 'off');
grid(ax2, 'on');
box(ax2, 'on');
xlabel(ax2, '時間 t（秒）');
ylabel(ax2, '追蹤誤差 e(t)');
title(ax2, sprintf('固定 H_r = %.4f（理想值），H_f 偏移 0%%、+1%%、+5%%、+10%%', Hr_nom));
legend(ax2, legendLabels(caseHf), 'Location', 'northeast');
styleAxes(ax2, cnFont);

% 穩態差異約 0.01，遠小於暫態峰值（約 0.25）。在主圖右下空白處放大末端。
zoomT0 = 20;
mask = t >= zoomT0;
yZoom = cell2mat(arrayfun(@(c) c.e(mask), caseHf, 'UniformOutput', false));
yMin = min(yZoom(:));
yMax = max(yZoom(:));
pad = max(0.25 * (yMax - yMin), 1e-3);
yLimZoom = [yMin - pad, yMax + pad];

rectangle(ax2, 'Position', [zoomT0, yLimZoom(1), tEnd - zoomT0, diff(yLimZoom)], ...
    'EdgeColor', [0.25, 0.25, 0.25], 'LineStyle', '--', 'LineWidth', 0.8);

% 依照主座標軸的實際位置，把放大圖放在圖二右下角的空白區
drawnow;
axPos = ax2.Position;
inset = axes('Parent', fig, 'Position', [ ...
    axPos(1) + axPos(3) * 0.56, ...
    axPos(2) + axPos(4) * 0.08, ...
    axPos(3) * 0.38, ...
    axPos(4) * 0.30]);
hold(inset, 'on');
for k = 1:nCase
    plot(inset, t, caseHf(k).e, 'Color', colors(k, :), ...
        'LineStyle', styles{k}, 'LineWidth', 1.4);
end
hold(inset, 'off');
xlim(inset, [zoomT0, tEnd]);
ylim(inset, yLimZoom);
grid(inset, 'on');
box(inset, 'on');
title(inset, sprintf('穩態放大（t = %g ~ %g s）', zoomT0, tEnd));
xlabel(inset, 't（秒）');
ylabel(inset, 'e(t)');
styleAxes(inset, cnFont);
set(inset, 'FontSize', 8);
inset.Title.FontSize = 9;

sg = sgtitle('參數公差對單位斜坡追蹤誤差的影響');
sg.FontName = cnFont;
sg.FontSize = 14;

%% 7. 儲存圖檔與數值，供後續報告使用
outPng = fullfile(scriptDir, 'type2_error_responses.png');
outFig = fullfile(scriptDir, 'type2_error_responses.fig');
try
    exportgraphics(fig, outPng, 'Resolution', 200);
catch
    saveas(fig, outPng);
end
savefig(fig, outFig);

save(fullfile(scriptDir, 'simulation_results.mat'), ...
    'A', 'tau', 'r_o', 'Hr_nom', 'Hf_nom', 'deltaPct', 't', 'r', ...
    'caseHr', 'caseHf');

fprintf('圖檔已寫入：%s\n', outPng);

%% ------------------------------------------------------------------------
function result = simulateError(A, tau, Hr, Hf, r, t, r_o, pct)
%SIMULATEERROR  由誤差轉移函數計算一組 (H_r, H_f) 的 e(t) 與漸近預測
%
% 分子：tau*s^2 + (1 - A*H_f)*s + A*(1 - H_r)
% 分母：tau*s^2 + s + A     （即 s*(tau*s + 1) + A）
% 輸入是斜坡而非步階，所以用 lsim，而不是直接對此轉移函數呼叫 step。

    num = [tau, 1 - A * Hf, A * (1 - Hr)];
    den = [tau, 1, A];
    sys = tf(num, den);

    e = lsim(sys, r, t);

    % 斜率項來自分子的常數項 A(1-H_r)；常數項來自分子的 s 項與分母常數項。
    slope  = (1 - Hr) * r_o;
    offset = r_o * (Hr / A - Hf);

    result.pct     = pct;
    result.Hr      = Hr;
    result.Hf      = Hf;
    result.e       = e;
    result.slope   = slope;
    result.offset  = offset;
    result.bounded = abs(slope) < 1e-10;
    result.ePred   = slope * t(end) + offset;
    result.eSim    = e(end);

    if result.bounded
        result.ess = offset;          % 等於 r_o * (1/A - H_f)
    elseif slope > 0
        result.ess = inf;
    else
        result.ess = -inf;
    end
end

function printScenario(fid, cases)
%PRINTSCENARIO  印出單一情境的參數、理論穩態特性與 t_end 數值

    logf(fid, '%6s  %10s  %10s  %-26s  %14s  %14s\n', ...
        '偏移', 'H_r', 'H_f', '理論穩態特性', '漸近式 e(30)', '模擬 e(30)');
    logf(fid, '%s\n', repmat('-', 1, 92));

    for k = 1:numel(cases)
        c = cases(k);
        if c.bounded
            trait = sprintf('有限  e_ss = % .6f', c.ess);
        else
            trait = sprintf('無界  斜率 = % .6f', c.slope);
        end
        logf(fid, '%5.0f%%  %10.6f  %10.6f  %-26s  %14.6f  %14.6f\n', ...
            c.pct, c.Hr, c.Hf, trait, c.ePred, c.eSim);
    end
end

function labels = legendLabels(cases)
%LEGENDLABELS  圖例同時標出偏移百分比與模擬終點值

    labels = cell(1, numel(cases));
    for k = 1:numel(cases)
        c = cases(k);
        if c.pct == 0
            tag = '理想 0%';
        else
            tag = sprintf('+%g%%', c.pct);
        end
        labels{k} = sprintf('%s   e(30) = %.4f', tag, c.eSim);
    end
end

function styleAxes(ax, fontName)
%STYLEAXES  統一座標軸中文字型，避免預設字型無法顯示中文

    set(ax, 'FontName', fontName, 'FontSize', 11, 'LineWidth', 0.8);
    ax.XLabel.FontName = fontName;
    ax.YLabel.FontName = fontName;
    ax.Title.FontName  = fontName;
    ax.Title.FontSize  = 12;
    if ~isempty(ax.Legend)
        ax.Legend.FontName = fontName;
        ax.Legend.FontSize = 9;
        ax.Legend.Interpreter = 'none';
        ax.Legend.AutoUpdate = 'off';
    end
end

function logf(fid, fmt, varargin)
%LOGF  同一段文字同時送到 Command Window 與結果檔

    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end
