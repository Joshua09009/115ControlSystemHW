%% HW05 Problem 3 — 直流馬達轉速前饋控制的延伸模擬（Stage 1）
% -------------------------------------------------------------------------
% 題目：Problem 4.38  Feedforward Control of a DC Motor Speed-Control System
%
% 本腳本只負責產生可重現的數值與圖，不撰寫最終討論。
% 請在 MATLAB 中執行：
%
%   cd('本檔所在的 hw05 資料夾')
%   hw05_p3_extension
%
% 模擬使用可控標準型狀態空間，並以矩陣指數做零階保持。
% 只需要 MATLAB 基底函數（expm、eig、roots），不需要 Control System Toolbox。
% 步階與取樣間隔內為常數的輸入，在取樣瞬間與連續時間響應一致。
%
% 標稱受控體與前饋：
%   G_nom(s) = b / (s^2 + a*s + c)
%   Kff      = G_nom(0)^(-1) = c / b
%   比例控制 kp = 3
%
% 控制律（r 為參考、w 為輸出端干擾、y 為轉速）：
%   純回授     u = Dc(s) * (r - y)
%   靜態前饋   u = Dc(s) * (r - y) + Kff * r - Kff * w
% 其中 Dc(s) = kp            （P）
%      Dc(s) = kp + ki/s     （PI）
%
% 輸出端干擾的定義是 y = w + G(s)*u。w 與 r 同時為 1 時，
% y = 1、u = 0 對任何穩定的 Dc 都是平衡點（因為 w 自己就補上了目標轉速）。
% Module 3 的 t = 5 s 誤差因此理論上全部接近 0，架構差異出現在暫態。
%
% 積分增益預設 ki = 4，而不是 15。kp 固定為 3 時，ki = 15 仍穩定，
% 但共軛極點阻尼比約 0.08，步階超越量很大，無法在 t = 2.5 s 之後的
% 2.5 秒視窗內以良好阻尼恢復。Module 3 會印出 ki = 2、4、15 的極點對照。
% 若要改回 15，只改下方參數區的 ki。
% -------------------------------------------------------------------------

clear; clc; close all;

%% 參數區（要做對照實驗時，優先改這裡）
b  = 59.292;     % 標稱分子
a  = 6.978;      % 標稱一次項
c  = 15.123;     % 標稱常數項，也是直流增益的分母
kp = 3;
ki = 4;

sigmaN  = 0.05;     % 干擾量測雜訊的標準差（每個取樣點獨立）
wc      = 20;       % 量測低通的截止頻率，rad/s
tEnd    = 5;        % 模擬終點，秒
dt      = 1e-3;     % 取樣間距，秒
tDisturb = 2.5;     % Module 3 干擾步階進入時刻
tProbe1  = 2.4;     % 追蹤階段讀數時刻（干擾進入前）
tProbe2  = 5.0;     % 干擾後讀數時刻
tJitter  = 3.0;     % 此時間之後視為穩態，用來計算抖動
ssBand   = 0.02;    % 安定／恢復時間的絕對誤差帶（單位步階的 2%）
rngSeed  = 20261022;

deltaList = [-0.20, -0.10, -0.05, 0, 0.05, 0.10, 0.20];
deltaShow = 0.15;   % Module 3 的增益偏差

G0  = b / c;        % 標稱直流增益
Kff = c / b;        % 標稱直流增益的倒數，精確使用 c/b 而不是四捨五入值

nStep = round(tEnd / dt);
t = (0:nStep).' * dt;

scriptDir = fileparts(mfilename('fullpath'));
figDir    = fullfile(scriptDir, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
cnFont = 'Microsoft JhengHei';

outTxt = fullfile(scriptDir, 'simulation_results.txt');
fid = fopen(outTxt, 'w', 'n', 'UTF-8');
if fid < 0
    error('無法建立結果檔：%s', outTxt);
end

%% 共用輸入
rStep = ones(size(t));
wStep = ones(size(t));
rZero = zeros(size(t));
wZero = zeros(size(t));

%% =====================================================================
banner(fid, '參數與標稱閉迴路');
logf(fid, 'G_nom(s) = %.6f / (s^2 + %.6f s + %.6f)\n', b, a, c);
logf(fid, 'G_nom(0) = %.8f\n', G0);
logf(fid, 'Kff = G_nom(0)^(-1) = %.8f   （四捨五入參考值 0.2550597）\n', Kff);
logf(fid, 'kp = %.6g，ki = %.6g，模擬 dt = %.4g s，tEnd = %.4g s\n', kp, ki, dt, tEnd);
logf(fid, '前饋相對純 P 的高頻／全頻縮放（同一受控體、純增益 Dc）：(kp+Kff)/kp = %.6f\n', (kp + Kff) / kp);
logf(fid, '追蹤穩態誤差（有前饋、增益偏差 δ）：e_ss = -δ / (1 + kp*(1+δ)*G_nom(0))\n');
logf(fid, '干擾穩態輸出（有前饋）：y_ss = e_ss。兩者共用同一個直流公式。\n');

AclNom = loopAcl(b, a, c, 0, kp, 0, false);
logf(fid, '\n標稱 P + 前饋的閉迴路極點（前饋不改變特徵根，P 與 P+前饋相同）：\n');
logPoles(fid, AclNom);

%% =====================================================================
banner(fid, 'Module 0  基準驗證：完美模型 G_act = G_nom');
logf(fid, '比較「純 P」與「P + 靜態前饋」。誤差帶 = ±%.3f（單位步階的 2%%）。\n', ssBand);
logf(fid, '安定時間的目標是各曲線自己的理論穩態值。\n');
logf(fid, '另外用「能否進入 |y-1| < %.3f」判斷純 P 是否真的跟上參考。\n', ssBand);

[yTrP,  uTrP]  = simCase(b, a, c, 0, kp, 0, false, 0,   0,   rStep, wZero, dt);
[yTrF,  uTrF]  = simCase(b, a, c, 0, kp, 0, false, Kff, Kff, rStep, wZero, dt);
[yDsP,  uDsP]  = simCase(b, a, c, 0, kp, 0, false, 0,   0,   rZero, wStep, dt);
[yDsF,  uDsF]  = simCase(b, a, c, 0, kp, 0, false, Kff, Kff, rZero, wStep, dt);

essP_thy = 1 / (1 + kp * G0);
essF_thy = 0;
yDsP_thy = essP_thy;
yDsF_thy = 0;

mTrP = responseMetrics(t, yTrP, 1 - essP_thy, ssBand);
mTrF = responseMetrics(t, yTrF, 1 - essF_thy, ssBand);
mDsP = responseMetrics(t, yDsP, yDsP_thy, ssBand);
mDsF = responseMetrics(t, yDsF, yDsF_thy, ssBand);

% 純 P 對「參考 = 1」的安定時間。穩態誤差大於誤差帶時，這個時間不存在。
TsRefP = settlingTime(t, yTrP, 1, ssBand);
TsRefF = settlingTime(t, yTrF, 1, ssBand);

scale = (kp + Kff) / kp;
scaleErr = max(abs(yTrF - scale * yTrP));

logf(fid, '\n--- 參考追蹤，r = 1(t)，w = 0 ---\n');
logf(fid, '%-16s  %12s  %12s\n', '項目', '純 P', 'P+前饋');
printPair(fid, '理論 e_ss',        essP_thy,          essF_thy,          '%12.6f');
printPair(fid, '模擬 e_ss=1-y(end)', 1 - yTrP(end),   1 - yTrF(end),     '%12.6e');
printPair(fid, '模擬 y(end)',      yTrP(end),         yTrF(end),         '%12.6f');
printPair(fid, '理論 y_ss',        1 - essP_thy,      1 - essF_thy,      '%12.6f');
printPair(fid, '峰值 y_peak',      mTrP.yPeak,        mTrF.yPeak,        '%12.6f');
printPair(fid, '峰值時刻 (s)',     mTrP.tPeak,        mTrF.tPeak,        '%12.4f');
printPair(fid, '超越量 相對穩態 (%)', mTrP.overshoot, mTrF.overshoot,    '%12.4f');
printPair(fid, '超過參考的量 (%)', mTrP.overRef,      mTrF.overRef,      '%12.4f');
printPair(fid, 'T_s 對自身穩態 (s)', mTrP.Ts,         mTrF.Ts,           '%12.4f');
printPair(fid, 'IAE=∫|1-y|dt',     trapz(t, abs(1 - yTrP)), trapz(t, abs(1 - yTrF)), '%12.6f');
printPair(fid, 'u(end)',           uTrP(end),         uTrF(end),         '%12.6f');
printPair(fid, 'max|u|',           max(abs(uTrP)),    max(abs(uTrF)),    '%12.6f');
logf(fid, '進入 |y-1|<%.3f 的時間 (s)：純 P = %s，P+前饋 = %s\n', ...
    ssBand, fmtMaybe(TsRefP, '%.4f'), fmtMaybe(TsRefF, '%.4f'));
logf(fid, '追蹤響應縮放檢核：max |y_ff - ((kp+Kff)/kp)*y_P| = %.3e\n', scaleErr);
logf(fid, '（此值只應有數值誤差。前饋不改變極點，只把純 P 的 y(t) 乘上常數。）\n');

logf(fid, '\n--- 輸出干擾，r = 0，w = 1(t) ---\n');
logf(fid, '%-16s  %12s  %12s\n', '項目', '純 P', 'P+前饋');
printPair(fid, '理論 y_ss',        yDsP_thy,          yDsF_thy,          '%12.6f');
printPair(fid, '模擬 y(end)',      yDsP(end),         yDsF(end),         '%12.6e');
printPair(fid, 'y(0+)',            yDsP(1),           yDsF(1),           '%12.6f');
printPair(fid, 'y_min',            mDsP.yMin,         mDsF.yMin,         '%12.6f');
printPair(fid, 'y_min 時刻 (s)',   mDsP.tMin,         mDsF.tMin,         '%12.4f');
printPair(fid, 'T_s 對理論穩態 (s)', mDsP.Ts,         mDsF.Ts,           '%12.4f');
printPair(fid, 'IAE=∫|y|dt',       trapz(t, abs(yDsP)), trapz(t, abs(yDsF)), '%12.6f');
printPair(fid, 'u(0+)',            uDsP(1),           uDsF(1),           '%12.6f');
printPair(fid, 'u(end)',           uDsP(end),         uDsF(end),         '%12.6f');
printPair(fid, 'max|u|',           max(abs(uDsP)),    max(abs(uDsF)),    '%12.6f');

fig0 = newFig('Module 0 基準驗證', 1080, 820, cnFont);
tl0 = tiledlayout(fig0, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl0);
plotPair(ax, t, yTrP, yTrF, '輸出 y(t)');
yline(ax, 1, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
title(ax, '參考追蹤');
xlim(ax, [0, 3]);
styleAxes(ax, cnFont);

ax = nexttile(tl0);
plotPair(ax, t, uTrP, uTrF, '控制輸入 u(t)');
title(ax, '追蹤時的控制力');
xlim(ax, [0, 3]);
styleAxes(ax, cnFont);

ax = nexttile(tl0);
plotPair(ax, t, yDsP, yDsF, '輸出 y(t)');
yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
title(ax, '干擾抑制');
xlim(ax, [0, 3]);
styleAxes(ax, cnFont);

ax = nexttile(tl0);
plotPair(ax, t, uDsP, uDsF, '控制輸入 u(t)');
title(ax, '抑制干擾時的控制力');
xlim(ax, [0, 3]);
styleAxes(ax, cnFont);
sg = title(tl0, 'Module 0  完美模型：純 P 與 P + 靜態前饋');
sg.FontName = cnFont;
sg.FontSize = 14;
sg.Interpreter = 'none';
saveFigure(fig0, figDir, 'fig0_baseline');

%% =====================================================================
banner(fid, 'Module 1  增益偏差：G_act = (1+δ) G_nom，Kff 固定為標稱倒數');
logf(fid, '前饋增益鎖在標稱值。δ 只放大真實馬達的分子，極點多項式不變。\n');
logf(fid, '超越量定義：(y_peak - y_ss理論) / |y_ss理論| × 100%%，且小於 0 時記為 0。\n');
logf(fid, '安定時間：|y(t) - y_ss理論| 此後都落在 %.3f 以內的第一個時刻。\n', ssBand);
logf(fid, '理論 e_ss = -δ / (1 + kp*(1+δ)*G_nom(0))。有前饋時，干擾的 y_ss 等於這個 e_ss。\n');

nDelta = numel(deltaList);
mod1(nDelta) = struct(); %#ok<NASGU>
for k = 1:nDelta
    delta = deltaList(k);
    [yT, uT] = simCase(b, a, c, delta, kp, 0, false, Kff, Kff, rStep, wZero, dt);
    [yW, uW] = simCase(b, a, c, delta, kp, 0, false, Kff, Kff, rZero, wStep, dt);

    essThy = -delta / (1 + kp * (1 + delta) * G0);
    ySsT   = 1 - essThy;
    ySsW   = essThy;
    metT   = responseMetrics(t, yT, ySsT, ssBand);
    metW   = responseMetrics(t, yW, ySsW, ssBand);

    mod1(k).delta   = delta;
    mod1(k).yT      = yT;
    mod1(k).uT      = uT;
    mod1(k).yW      = yW;
    mod1(k).uW      = uW;
    mod1(k).essThy  = essThy;
    mod1(k).essSim  = 1 - yT(end);
    mod1(k).os      = metT.overshoot;
    mod1(k).TsT     = metT.Ts;
    mod1(k).yPeak   = metT.yPeak;
    mod1(k).ySsT    = yT(end);
    mod1(k).ySsWThy = ySsW;
    mod1(k).ySsW    = yW(end);
    mod1(k).yMinW   = metW.yMin;
    mod1(k).TsW     = metW.Ts;
end

logf(fid, '\n--- 追蹤（有前饋） ---\n');
logf(fid, '%7s  %12s  %12s  %12s  %12s  %10s  %10s\n', ...
    'δ(%)', 'e_ss理論', 'e_ss模擬', '理論-模擬', 'y_peak', 'OS(%)', 'T_s(s)');
logf(fid, '%s\n', repmat('-', 1, 86));
for k = 1:nDelta
    logf(fid, '%+7.1f  %+12.6f  %+12.6e  %+12.3e  %12.6f  %10.4f  %10s\n', ...
        100 * mod1(k).delta, mod1(k).essThy, mod1(k).essSim, ...
        mod1(k).essThy - mod1(k).essSim, mod1(k).yPeak, mod1(k).os, ...
        fmtMaybe(mod1(k).TsT, '%.4f'));
end

logf(fid, '\n--- 干擾抑制（有前饋）。理論 y_ss 應等於上一表的 e_ss理論 ---\n');
logf(fid, '%7s  %12s  %12s  %12s  %12s  %10s\n', ...
    'δ(%)', 'y_ss理論', 'y_ss模擬', '理論-模擬', 'y_min', 'T_s(s)');
logf(fid, '%s\n', repmat('-', 1, 74));
for k = 1:nDelta
    logf(fid, '%+7.1f  %+12.6f  %+12.6e  %+12.3e  %12.6f  %10s\n', ...
        100 * mod1(k).delta, mod1(k).ySsWThy, mod1(k).ySsW, ...
        mod1(k).ySsWThy - mod1(k).ySsW, mod1(k).yMinW, ...
        fmtMaybe(mod1(k).TsW, '%.4f'));
end

essErrMax = max(abs([mod1.essThy] - [mod1.essSim]));
yWErrMax  = max(abs([mod1.ySsWThy] - [mod1.ySsW]));
logf(fid, '\n追蹤 |e_ss理論 - e_ss模擬| 的最大值 = %.3e\n', essErrMax);
logf(fid, '干擾 |y_ss理論 - y_ss模擬| 的最大值 = %.3e\n', yWErrMax);

fig1 = newFig('Module 1 增益偏差', 1080, 820, cnFont);
tl1 = tiledlayout(fig1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
colors1 = mismatchColors(nDelta);
labels1 = cell(1, nDelta);
for k = 1:nDelta
    labels1{k} = sprintf('δ = %+.0f%%', 100 * deltaList(k));
end

ax = nexttile(tl1);
hold(ax, 'on');
for k = 1:nDelta
    lw = 1.3 + 0.9 * (deltaList(k) == 0);
    plot(ax, t, mod1(k).yT, 'Color', colors1(k, :), 'LineWidth', lw);
end
yline(ax, 1, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, 3]);
xlabel(ax, '時間 t（秒）');
ylabel(ax, '輸出 y(t)');
title(ax, '步階參考追蹤');
legend(ax, labels1, 'Location', 'eastoutside');
styleAxes(ax, cnFont);

ax = nexttile(tl1);
hold(ax, 'on');
for k = 1:nDelta
    lw = 1.3 + 0.9 * (deltaList(k) == 0);
    plot(ax, t, mod1(k).yW, 'Color', colors1(k, :), 'LineWidth', lw);
end
yline(ax, 0, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, 3]);
xlabel(ax, '時間 t（秒）');
ylabel(ax, '輸出 y(t)');
title(ax, '步階干擾抑制');
legend(ax, labels1, 'Location', 'eastoutside');
styleAxes(ax, cnFont);
sg = title(tl1, 'Module 1  Kff 固定，真實增益 Gact = (1+δ) Gnom');
sg.FontName = cnFont;
sg.FontSize = 14;
sg.Interpreter = 'none';
saveFigure(fig1, figDir, 'fig1_gain_mismatch');

%% =====================================================================
banner(fid, 'Module 2  干擾初值突波、負向超越，以及量測雜訊');
logf(fid, '受控體為標稱模型。干擾是加在輸出端的單位步階，不是加在 u 上。\n');
logf(fid, 'G(s) 的相對階數 = 2，所以任何有限的 u 都無法在 t = 0+ 改變 y。\n');

% 完美模型下，干擾轉移函數分子 = s^2 + a s + (c - Kff*b) = s(s+a)
logf(fid, '\n--- 初值（有前饋，標稱模型） ---\n');
logf(fid, '理論 y(0+) = lim s→∞ T_w(s) = 1\n');
logf(fid, '模擬 y(0)          = %.6f\n', yDsF(1));
logf(fid, '理論 y''(0+) = 0     （步階干擾本身斜率為 0，受控體相對階數 2，u 尚無貢獻）\n');
u0Thy = -(kp + Kff);
ydd0  = b * u0Thy;
logf(fid, '理論 u(0+) = -(kp + Kff) = %.6f\n', u0Thy);
logf(fid, '模擬 u(0)          = %.6f\n', uDsF(1));
logf(fid, '理論 y''''(0+) = b*u(0+) = %.6f   （離開 1 的初始曲率為負）\n', ydd0);
logf(fid, '理論 y(∞) = 0，理論 u(∞) = -Kff = %.6f\n', -Kff);
logf(fid, '模擬 y_min = %.6f，發生在 t = %.4f s\n', mDsF.yMin, mDsF.tMin);
logf(fid, '模擬 u(end) = %.6f，y(end) = %.3e\n', uDsF(end), yDsF(end));
logf(fid, '無前饋對照：y(0+) = %.6f（同樣是 1），y_min = %.6f，u(0+) = %.6f，y(∞) = %.6f\n', ...
    yDsP(1), mDsP.yMin, uDsP(1), yDsP(end));
logf(fid, '標稱閉迴路極點與阻尼：\n');
logPoles(fid, AclNom);
zTw = roots([1, a, c - Kff * b]);
logf(fid, 'T_w(s) 的零點（分子 s^2 + a s + (c - Kff*b)）：\n');
for iz = 1:numel(zTw)
    logf(fid, '    z = %+.6f%+.6fj\n', real(zTw(iz)), imag(zTw(iz)));
end

fig2a = newFig('Module 2 突波', 1080, 720, cnFont);
tl2a = tiledlayout(fig2a, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl2a);
plotPair(ax, t, yDsP, yDsF, '輸出 y(t)');
yline(ax, 0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
yline(ax, 1, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
hold(ax, 'on');
plot(ax, mDsF.tMin, mDsF.yMin, 'o', 'Color', [0.00 0.35 0.70], ...
    'MarkerSize', 7, 'LineWidth', 1.4, 'HandleVisibility', 'off');
hold(ax, 'off');
title(ax, sprintf('輸出突波：y(0+) = 1，有前饋時 y_min = %.3f', mDsF.yMin));
xlim(ax, [0, 3]);
styleAxes(ax, cnFont);

ax = nexttile(tl2a);
plotPair(ax, t, uDsP, uDsF, '控制輸入 u(t)');
title(ax, sprintf('控制力：有前饋時 u(0+) = %.3f，u(∞) = %.3f', uDsF(1), -Kff));
xlim(ax, [0, 3]);
styleAxes(ax, cnFont);
sg = title(tl2a, 'Module 2  輸出端步階干擾的初值與控制力（標稱模型）');
sg.FontName = cnFont;
sg.FontSize = 14;
sg.Interpreter = 'none';
saveFigure(fig2a, figDir, 'fig2_disturbance_spike');

% ----- 量測雜訊 -----
% 控制器看到的是 w_hat，真實干擾仍是 w。回授使用真實 y，雜訊只進前饋路徑。
logf(fid, '\n--- 干擾量測雜訊 ---\n');
logf(fid, 'w_hat = w + n，n 為每個取樣點獨立的 N(0, σ^2)，σ = %.3f，rng = %d。\n', sigmaN, rngSeed);
logf(fid, '這是取樣週期 %.4g s 上的離散白雜訊，不是單位頻譜密度的連續白雜訊。\n', dt);
logf(fid, '低通 H(s) = wc/(s+wc)，wc = %.4g rad/s（%.3f Hz），時間常數 %.4f s。\n', ...
    wc, wc / (2 * pi), 1 / wc);
logf(fid, '高頻時 u/n → -Kff / (1 + kp*G) ≈ -Kff = %.6f，雜訊被前饋增益縮小，不是被放大。\n', Kff);
logf(fid, '低頻時 |u/n| 再被 |1+kp*G(0)| = %.4f 縮小。\n', abs(1 + kp * G0));
logf(fid, '全時域 RMS(u) 會被 t = 0 的大控制力主導；請同時看 t ≥ %.1f s 的穩態抖動 RMS。\n', tJitter);

rng(rngSeed, 'twister');
n = sigmaN * randn(size(t));
wHatRaw = wStep + n;
wHatLp  = lowpassZOH(wc, wHatRaw, dt);

[yClean, uClean] = simMeasured(b, a, c, kp, Kff, wStep, wStep,   dt);
[yRaw,   uRaw]   = simMeasured(b, a, c, kp, Kff, wStep, wHatRaw, dt);
[yLp,    uLp]    = simMeasured(b, a, c, kp, Kff, wStep, wHatLp,  dt);
logf(fid, '自檢：無雜訊前饋與 Module 0 干擾響應的最大差 |y| = %.3e，|u| = %.3e\n', ...
    max(abs(yClean - yDsF)), max(abs(uClean - uDsF)));

noiseCases = struct( ...
    'name', {'無雜訊', '雜訊、無濾波', '雜訊 + 低通'}, ...
    'y',    {yClean, yRaw, yLp}, ...
    'u',    {uClean, uRaw, uLp});

logf(fid, '\n%16s  %12s  %12s  %14s  %14s  %12s\n', ...
    '案例', 'RMS(u)全時域', 'max|u|', '穩態抖動RMS(u)', '穩態抖動RMS(y)', 'y_min');
logf(fid, '%s\n', repmat('-', 1, 86));
for k = 1:numel(noiseCases)
    uk = noiseCases(k).u;
    yk = noiseCases(k).y;
    logf(fid, '%16s  %12.6f  %12.6f  %14.6f  %14.6f  %12.6f\n', ...
        noiseCases(k).name, rmsValue(uk), max(abs(uk)), ...
        jitterRms(t, uk, tJitter), jitterRms(t, yk, tJitter), min(yk));
end
logf(fid, '未濾波時，穩態 u 抖動的粗估 Kff*σ = %.6f（實際值還會再被靈敏度函數濾過）。\n', Kff * sigmaN);

fig2b = newFig('Module 2 雜訊', 1080, 820, cnFont);
tl2b = tiledlayout(fig2b, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
noiseColor = [ ...
    0.15 0.15 0.15
    0.85 0.33 0.10
    0.00 0.45 0.70];
noiseStyle = {'-', '-', '-'};
noiseWidth = [2.0, 1.0, 1.6];
noiseLabel = {'無雜訊', '雜訊、無濾波', '雜訊 + 低通'};

ax = nexttile(tl2b);
plotCases(ax, t, {yClean, yRaw, yLp}, noiseColor, noiseStyle, noiseWidth, noiseLabel);
ylabel(ax, '輸出 y(t)');
title(ax, '全時域輸出');
xlim(ax, [0, tEnd]);
styleAxes(ax, cnFont);

ax = nexttile(tl2b);
plotCases(ax, t, {uClean, uRaw, uLp}, noiseColor, noiseStyle, noiseWidth, noiseLabel);
ylabel(ax, '控制輸入 u(t)');
title(ax, '全時域控制力');
xlim(ax, [0, tEnd]);
styleAxes(ax, cnFont);

ax = nexttile(tl2b);
plotZoom(ax, t, [yClean, yRaw, yLp], tJitter, tEnd, noiseColor, noiseStyle, noiseWidth);
ylabel(ax, '輸出 y(t)');
title(ax, sprintf('穩態放大  t = %.0f–%.0f s', tJitter, tEnd));
styleAxes(ax, cnFont);

ax = nexttile(tl2b);
plotZoom(ax, t, [uClean, uRaw, uLp], tJitter, tEnd, noiseColor, noiseStyle, noiseWidth);
ylabel(ax, '控制輸入 u(t)');
title(ax, sprintf('穩態放大  t = %.0f–%.0f s', tJitter, tEnd));
styleAxes(ax, cnFont);
sg = title(tl2b, sprintf('Module 2  干擾量測雜訊 σ = %.2f，低通 wc = %.0f rad/s', sigmaN, wc));
sg.FontName = cnFont;
sg.FontSize = 14;
sg.Interpreter = 'none';
saveFigure(fig2b, figDir, 'fig2_sensor_noise');

%% =====================================================================
banner(fid, 'Module 3  四種架構，δ = +15%，先追蹤再加干擾');
logf(fid, '真實受控體 G_act = (1+%.2f) G_nom。Kff 仍使用標稱倒數。\n', deltaShow);
logf(fid, 'r(t) = 1，t ≥ 0；w(t) = 1，t ≥ %.2f s。圖的時間軸為 [0, %.0f] s。\n', tDisturb, tEnd);
logf(fid, '讀數時刻：追蹤階段 t = %.2f s（實際取樣 t = %.4f s），干擾後 t = %.2f s。\n', ...
    tProbe1, t(round(tProbe1 / dt) + 1), tProbe2);
logf(fid, ['前饋架構的控制律是 (a)(b) 合併：u = Dc(r-y) + Kff*r - Kff*w。\n' ...
    '純 P、純 PI 不含這兩項前饋。\n']);
logf(fid, ['當 r(∞) = w(∞) = 1 時，y = 1 且 u = 0 是所有四種架構的平衡點，\n' ...
    '所以 t = %.1f s 的誤差理論上都是 0。比較重點是 t = %.1f s 的追蹤誤差、\n' ...
    '追蹤超越量，以及干擾之後回到 |y-1| < %.3f 的恢復時間。\n'], ...
    tProbe2, tProbe1, ssBand);

logf(fid, '\n--- ki 候選（標稱受控體、純 PI，kp = %.4g） ---\n', kp);
logf(fid, '特徵式 s^3 + a s^2 + (c + b*kp) s + b*ki。Routh 條件 a*(c+b*kp) > b*ki。\n');
logf(fid, '%8s  %14s  %10s  %10s\n', 'ki', 'Routh裕度', 'zeta', '實極點');
kiCand = [2, 4, 15];
for k = 1:numel(kiCand)
    info = piDesignInfo(a, b, c, kp, kiCand(k));
    mark = '';
    if abs(kiCand(k) - ki) < 1e-12
        mark = '  <-- 本腳本使用';
    end
    logf(fid, '%8.2f  %14.4f  %10.4f  %10.4f%s\n', ...
        kiCand(k), info.routh, info.zeta, info.realPole, mark);
end
logf(fid, '本腳本選定 ki = %.4g。ki = 15 穩定但阻尼明顯低於純 P（純 P 的 zeta 約 0.25）。\n', ki);

Gmis0 = (1 + deltaShow) * G0;
essP_mis  = 1 / (1 + kp * Gmis0);
essFF_mis = (1 - Kff * Gmis0) / (1 + kp * Gmis0);
uSsP_mis  = kp * essP_mis;
uSsFF_mis = kp * essFF_mis + Kff;
uSsPI_mis = 1 / Gmis0;

arch(1).name = '純 P';
arch(1).useInt = false;
arch(1).Kr   = 0;
arch(1).Kw   = 0;
arch(1).e1Thy = essP_mis;
arch(1).u1Thy = uSsP_mis;

arch(2).name = 'P+靜態前饋';
arch(2).useInt = false;
arch(2).Kr   = Kff;
arch(2).Kw   = Kff;
arch(2).e1Thy = essFF_mis;
arch(2).u1Thy = uSsFF_mis;

arch(3).name = '純 PI';
arch(3).useInt = true;
arch(3).Kr   = 0;
arch(3).Kw   = 0;
arch(3).e1Thy = 0;
arch(3).u1Thy = uSsPI_mis;

arch(4).name = '二自由度 PI+前饋';
arch(4).useInt = true;
arch(4).Kr   = Kff;
arch(4).Kw   = Kff;
arch(4).e1Thy = 0;
arch(4).u1Thy = uSsPI_mis;

wShow = zeros(size(t));
wShow(round(tDisturb / dt) + 1:end) = 1;
rShow = rStep;

logf(fid, '\n相位一（w = 0）理論：e_ss(純 P) = %.6f，e_ss(P+前饋) = %.6f，e_ss(PI) = 0。\n', ...
    essP_mis, essFF_mis);
logf(fid, '相位一理論 u_ss：純 P = %.6f，P+前饋 = %.6f，PI 與 PI+前饋 = 1/G_act(0) = %.6f。\n', ...
    uSsP_mis, uSsFF_mis, uSsPI_mis);
logf(fid, '相位二（r = w = 1）理論：四種架構皆為 e_ss = 0、u_ss = 0。\n');

logf(fid, '\n失配受控體上的閉迴路極點（前饋不改變極點，所以只有兩組）：\n');
logf(fid, '純 P / P+前饋：\n');
logPoles(fid, loopAcl(b, a, c, deltaShow, kp, 0, false));
logf(fid, '純 PI / PI+前饋：\n');
logPoles(fid, loopAcl(b, a, c, deltaShow, kp, ki, true));
logf(fid, '追蹤通道 Y/R 的有限零點：純 PI 在 s = %.4f，PI+前饋在 s = %.4f。\n', ...
    -ki / kp, -ki / (kp + Kff));
logf(fid, '（PI 零點 = -ki/kp；PI+前饋零點 = -ki/(kp+Kff)。前饋只移動零點。）\n');

pre = t < tDisturb;
for k = 1:numel(arch)
    kiUse = 0;
    if arch(k).useInt
        kiUse = ki;
    end
    [yk, uk] = simCase(b, a, c, deltaShow, kp, kiUse, arch(k).useInt, ...
        arch(k).Kr, arch(k).Kw, rShow, wShow, dt);
    [e1, t1] = sampleAt(t, 1 - yk, tProbe1, dt);
    [e2, t2] = sampleAt(t, 1 - yk, tProbe2, dt);
    [u1, ~]  = sampleAt(t, uk, tProbe1, dt);
    [u2, ~]  = sampleAt(t, uk, tProbe2, dt);
    yPre = yk(pre);
    yPeak = max(yPre);
    ySs1 = 1 - arch(k).e1Thy;
    if abs(ySs1) > 1e-8
        os = max(0, (yPeak - ySs1) / abs(ySs1) * 100);
    else
        os = NaN;
    end
    idxBefore = find(pre, 1, 'last');
    idxAfter  = idxBefore + 1;
    Tr = recoveryTime(t, yk, tDisturb, 1, ssBand);

    arch(k).y = yk;
    arch(k).u = uk;
    arch(k).e1 = e1;
    arch(k).e2 = e2;
    arch(k).t1 = t1;
    arch(k).t2 = t2;
    arch(k).u1 = u1;
    arch(k).u2 = u2;
    arch(k).yPeak = yPeak;
    arch(k).os = os;
    arch(k).yBefore = yk(idxBefore);
    arch(k).yAfter  = yk(idxAfter);
    arch(k).yMinPost = min(yk(~pre));
    arch(k).Tr = Tr;
    arch(k).maxAbsU = max(abs(uk));
end

logf(fid, '\n--- 追蹤階段（干擾進入前） ---\n');
logf(fid, '%18s  %12s  %12s  %12s  %12s  %10s  %12s\n', ...
    '架構', 'e_ss理論', 'e(2.4)', 'u_ss理論', 'u(2.4)', 'OS(%)', 'y_peak');
logf(fid, '%s\n', repmat('-', 1, 96));
for k = 1:numel(arch)
    logf(fid, '%18s  %+12.6f  %+12.6f  %+12.6f  %+12.6f  %10.4f  %12.6f\n', ...
        arch(k).name, arch(k).e1Thy, arch(k).e1, arch(k).u1Thy, arch(k).u1, ...
        arch(k).os, arch(k).yPeak);
end

logf(fid, '\n--- 干擾進入與恢復。恢復時間對準 |y-1| < %.3f，從 t = %.2f s 起算 ---\n', ...
    ssBand, tDisturb);
logf(fid, '%18s  %12s  %12s  %12s  %12s  %12s  %12s\n', ...
    '架構', 'y(2.5-)', 'y(2.5+)', 'e(5.0)', 'u(5.0)', 'y_min恢復', 'T_rec(s)');
logf(fid, '%s\n', repmat('-', 1, 100));
for k = 1:numel(arch)
    logf(fid, '%18s  %12.6f  %12.6f  %+12.6f  %+12.6f  %12.6f  %12s\n', ...
        arch(k).name, arch(k).yBefore, arch(k).yAfter, arch(k).e2, arch(k).u2, ...
        arch(k).yMinPost, fmtMaybe(arch(k).Tr, '%.4f'));
end

logf(fid, '\n--- 控制力峰值（全時域） ---\n');
for k = 1:numel(arch)
    logf(fid, '%18s  max|u| = %.6f\n', arch(k).name, arch(k).maxAbsU);
end

archColor = [ ...
    0.35 0.35 0.35
    0.00 0.35 0.70
    0.85 0.33 0.10
    0.00 0.50 0.25];
archStyle = {'--', '-', '--', '-'};
archWidth = [1.6, 1.8, 1.6, 1.8];
archLabel = {arch.name};

fig3 = newFig('Module 3 架構比較', 1080, 820, cnFont);
tl3 = tiledlayout(fig3, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl3);
hold(ax, 'on');
for k = 1:numel(arch)
    plot(ax, t, arch(k).y, 'Color', archColor(k, :), ...
        'LineStyle', archStyle{k}, 'LineWidth', archWidth(k));
end
xl = xline(ax, tDisturb, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0);
xl.HandleVisibility = 'off';
yline(ax, 1, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, tEnd]);
xlabel(ax, '時間 t（秒）');
ylabel(ax, '輸出 y(t)');
title(ax, sprintf('輸出。垂直虛線是 t = %.1f s 的單位干擾', tDisturb));
legend(ax, archLabel, 'Location', 'eastoutside');
styleAxes(ax, cnFont);

ax = nexttile(tl3);
hold(ax, 'on');
for k = 1:numel(arch)
    plot(ax, t, arch(k).u, 'Color', archColor(k, :), ...
        'LineStyle', archStyle{k}, 'LineWidth', archWidth(k));
end
xl = xline(ax, tDisturb, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.0);
xl.HandleVisibility = 'off';
yline(ax, 0, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
hold(ax, 'off');
grid(ax, 'on'); box(ax, 'on');
xlim(ax, [0, tEnd]);
xlabel(ax, '時間 t（秒）');
ylabel(ax, '控制輸入 u(t)');
title(ax, '控制力');
legend(ax, archLabel, 'Location', 'eastoutside');
styleAxes(ax, cnFont);
sg = title(tl3, sprintf('Module 3  δ = +%.0f%%，單位參考於 t = 0，單位干擾於 t = %.1f s', ...
    100 * deltaShow, tDisturb));
sg.FontName = cnFont;
sg.FontSize = 14;
sg.Interpreter = 'none';
saveFigure(fig3, figDir, 'fig3_architecture_showdown');

%% 存檔
save(fullfile(scriptDir, 'simulation_results.mat'), ...
    'b', 'a', 'c', 'kp', 'ki', 'Kff', 'G0', 'sigmaN', 'wc', ...
    't', 'dt', 'deltaList', 'deltaShow', 'ssBand', ...
    'yTrP', 'yTrF', 'uTrP', 'uTrF', 'yDsP', 'yDsF', 'uDsP', 'uDsF', ...
    'mod1', 'noiseCases', 'n', 'arch');

logf(fid, '\n圖檔與數值已寫入：\n');
logf(fid, '  %s\n', fullfile(figDir, 'fig0_baseline.png'));
logf(fid, '  %s\n', fullfile(figDir, 'fig1_gain_mismatch.png'));
logf(fid, '  %s\n', fullfile(figDir, 'fig2_disturbance_spike.png'));
logf(fid, '  %s\n', fullfile(figDir, 'fig2_sensor_noise.png'));
logf(fid, '  %s\n', fullfile(figDir, 'fig3_architecture_showdown.png'));
logf(fid, '  %s\n', outTxt);

fclose(fid);
fprintf('結果表已寫入：%s\n', outTxt);

%% =====================================================================
function [y, u] = simCase(b, aCoeff, c, delta, kp, ki, useInt, KffR, KffW, r, w, dt)
%SIMCASE  模擬 y = w + G*u，u = Dc(r-y) + KffR*r - KffW*w。
% useInt = false 時 Dc = kp；useInt = true 時 Dc = kp + ki/s。
    b1 = (1 + delta) * b;
    [Acl, Bcl, Cy, Dy, Cu, Du] = loopMatrices(b1, aCoeff, c, kp, ki, useInt, KffR, KffW);
    [y, u] = zohSim(Acl, Bcl, Cy, Dy, Cu, Du, [r(:), w(:)], dt);
end

function [y, u] = simMeasured(b, aCoeff, c, kp, Kff, w, wHat, dt)
%SIMMEASURED  r = 0，真實干擾為 w，前饋使用的量測為 wHat。
%   u = -kp*y - Kff*wHat，y = w + G*u。
    b1 = b;
    A = [0, 1; -c, -aCoeff];
    B = [0; b1];
    C = [1, 0];
    Acl = A - kp * (B * C);
    Bcl = [B * (-kp), B * (-Kff)];
    Cy = C;
    Dy = [1, 0];
    Cu = -kp * C;
    Du = [-kp, -Kff];
    [y, u] = zohSim(Acl, Bcl, Cy, Dy, Cu, Du, [w(:), wHat(:)], dt);
end

function Acl = loopAcl(b, aCoeff, c, delta, kp, ki, useInt)
%LOOPACL  閉迴路系統矩陣。靜態前饋不進入特徵根。
    [Acl, ~, ~, ~, ~, ~] = loopMatrices((1 + delta) * b, aCoeff, c, kp, ki, useInt, 0, 0);
end

function [Acl, Bcl, Cy, Dy, Cu, Du] = loopMatrices(b1, aCoeff, c, kp, ki, useInt, KffR, KffW)
%LOOPMATRICES  受控體可控標準型，輸入順序為 [r, w]。
%
%   x' = [0 1; -c -a] x + [0; b1] u
%   y  = x1 + w
%   純 P：u = kp*(r-y) + KffR*r - KffW*w
%   PI ：再加 ki*xi，xi' = r - y
    A = [0, 1; -c, -aCoeff];
    B = [0; b1];
    C = [1, 0];
    Ap = A - kp * (B * C);

    if ~useInt
        Acl = Ap;
        Bcl = [B * (kp + KffR), B * (-kp - KffW)];
        Cy  = C;
        Dy  = [0, 1];
        Cu  = -kp * C;
        Du  = [kp + KffR, -kp - KffW];
        return
    end

    Acl = [Ap, ki * B; -C, 0];
    Bcl = [B * (kp + KffR), B * (-kp - KffW); 1, -1];
    Cy  = [C, 0];
    Dy  = [0, 1];
    Cu  = [-kp * C, ki];
    Du  = [kp + KffR, -kp - KffW];
end

function [y, u] = zohSim(Acl, Bcl, Cy, Dy, Cu, Du, U, dt)
%ZOHSIM  零階保持下的精確取樣。U 的第 k 列在 [t_k, t_{k+1}) 保持常數。
% 輸出先用當下狀態與當下輸入計算，因此直接項（例如 y 對 w）出現在步階的同一點。
    n = size(Acl, 1);
    nu = size(Bcl, 2);
    M = expm([Acl, Bcl; zeros(nu, n + nu)] * dt);
    Ad = M(1:n, 1:n);
    Bd = M(1:n, n + 1:n + nu);

    N = size(U, 1);
    y = zeros(N, 1);
    u = zeros(N, 1);
    z = zeros(n, 1);
    for k = 1:N
        uk = U(k, :).';
        y(k) = Cy * z + Dy * uk;
        u(k) = Cu * z + Du * uk;
        if k < N
            z = Ad * z + Bd * uk;
        end
    end
end

function v = lowpassZOH(wc, raw, dt)
%LOWPASSZOH  因果一階低通 v' = -wc*v + wc*raw，v(0-) = 0，沒有直接項。
    [v, ~] = zohSim(-wc, wc, 1, 0, 1, 0, raw(:), dt);
end

function met = responseMetrics(t, y, yFinal, band)
%RESPONSEMETRICS  峰值、相對穩態值的超越量、最低點與安定時間。
    [yPeak, iPeak] = max(y);
    [yMin, iMin] = min(y);
    if abs(yFinal) > 1e-8
        overshoot = max(0, (yPeak - yFinal) / abs(yFinal) * 100);
    else
        overshoot = NaN;
    end
    met.yPeak = yPeak;
    met.tPeak = t(iPeak);
    met.yMin = yMin;
    met.tMin = t(iMin);
    met.overshoot = overshoot;
    met.overRef = max(0, (yPeak - 1) * 100);
    met.Ts = settlingTime(t, y, yFinal, band);
end

function Ts = settlingTime(t, y, yFinal, band)
%SETTLINGTIME  最後一次離開誤差帶的下一個取樣時刻。到終點仍在帶外則回 NaN。
    outside = abs(y(:) - yFinal) > band;
    idx = find(outside, 1, 'last');
    if isempty(idx)
        Ts = 0;
    elseif idx >= numel(t)
        Ts = NaN;
    else
        Ts = t(idx + 1);
    end
end

function Tr = recoveryTime(t, y, t0, yFinal, band)
%RECOVERYTIME  從 t0 起算，此後都留在誤差帶內所需的時間。
    mask = t >= t0;
    tt = t(mask);
    yy = y(mask);
    outside = abs(yy - yFinal) > band;
    idx = find(outside, 1, 'last');
    if isempty(idx)
        Tr = 0;
    elseif idx >= numel(tt)
        Tr = NaN;
    else
        Tr = tt(idx + 1) - t0;
    end
end

function [val, tUsed] = sampleAt(t, x, tq, dt)
%SAMPLEAT  取最接近 tq 的格點。時間軸是 (k)*dt，k 從 0 開始。
    idx = round(tq / dt) + 1;
    idx = max(1, min(numel(t), idx));
    val = x(idx);
    tUsed = t(idx);
end

function info = piDesignInfo(aCoeff, b, c, kp, ki)
%PIDESIGNINFO  標稱受控體上純 PI 的 Routh 裕度、共軛阻尼比與實極點。
    beta = c + b * kp;
    gamma = b * ki;
    info.routh = aCoeff * beta - gamma;
    p = roots([1, aCoeff, beta, gamma]);
    realPole = NaN;
    zeta = NaN;
    for i = 1:numel(p)
        if abs(imag(p(i))) < 1e-8
            realPole = real(p(i));
        else
            wn = abs(p(i));
            zeta = -real(p(i)) / wn;
        end
    end
    info.realPole = realPole;
    info.zeta = zeta;
end

function r = rmsValue(x)
    x = x(:);
    r = sqrt(mean(x .^ 2));
end

function r = jitterRms(t, x, t0)
%JITTERRMS  穩態區間內去掉平均值之後的 RMS，用來代表抖動而不是直流成分。
    x = x(t >= t0);
    x = x(:) - mean(x);
    r = sqrt(mean(x .^ 2));
end

function s = fmtMaybe(x, fmt)
    if isnan(x)
        s = '未安定';
    else
        s = sprintf(fmt, x);
    end
end

function logPoles(fid, Acl)
    p = eig(Acl);
    p = p(:);
    [~, ix] = sortrows([real(p), imag(p)]);
    p = p(ix);
    if isempty(p)
        logf(fid, '    （無極點）\n');
        return
    end
    for i = 1:numel(p)
        wn = abs(p(i));
        if wn < 1e-10
            zeta = NaN;
        else
            zeta = -real(p(i)) / wn;
        end
        stable = '左半平面';
        if real(p(i)) > 1e-8
            stable = '右半平面';
        elseif abs(real(p(i))) <= 1e-8
            stable = '虛軸';
        end
        logf(fid, '    p = %+9.4f %+9.4fj    wn = %7.4f    zeta = %7.4f    %s\n', ...
            real(p(i)), imag(p(i)), wn, zeta, stable);
    end
end

function colors = mismatchColors(n)
%MISMATCHCOLORS  負偏差偏藍、零偏差為黑、正偏差偏紅。
    anchors = [ ...
        0.12 0.30 0.65
        0.25 0.50 0.78
        0.55 0.72 0.86
        0.15 0.15 0.15
        0.90 0.62 0.35
        0.82 0.32 0.18
        0.65 0.12 0.12];
    if n == size(anchors, 1)
        colors = anchors;
        return
    end
    xq = linspace(1, size(anchors, 1), n);
    colors = interp1(1:size(anchors, 1), anchors, xq);
end

function fig = newFig(figName, w, h, fontName)
    fig = figure('Name', figName, 'Color', 'w', 'Position', [60, 30, w, h]);
    set(fig, 'DefaultAxesFontName', fontName, ...
             'DefaultTextFontName', fontName, ...
             'DefaultAxesFontSize', 11, ...
             'InvertHardcopy', 'off');
end

function plotPair(ax, t, yP, yF, yLabelText)
    hold(ax, 'on');
    plot(ax, t, yP, '--', 'Color', [0.40 0.40 0.40], 'LineWidth', 1.6);
    plot(ax, t, yF, '-',  'Color', [0.00 0.35 0.70], 'LineWidth', 1.8);
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '時間 t（秒）');
    ylabel(ax, yLabelText);
    legend(ax, {'純 P', 'P + 前饋'}, 'Location', 'best');
end

function plotCases(ax, t, series, colors, styles, widths, labels)
    hold(ax, 'on');
    for k = 1:numel(series)
        plot(ax, t, series{k}, 'Color', colors(k, :), ...
            'LineStyle', styles{k}, 'LineWidth', widths(k));
    end
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '時間 t（秒）');
    legend(ax, labels, 'Location', 'best');
end

function plotZoom(ax, t, Y, t0, t1, colors, styles, widths)
%PLOTZOOM  只畫穩態區間，縱軸才會放大到雜訊看得到的尺度。
    mask = (t >= t0) & (t <= t1);
    tt = t(mask);
    hold(ax, 'on');
    for k = 1:size(Y, 2)
        plot(ax, tt, Y(mask, k), 'Color', colors(k, :), ...
            'LineStyle', styles{k}, 'LineWidth', widths(k));
    end
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [t0, t1]);
    xlabel(ax, '時間 t（秒）');
    yy = Y(mask, :);
    span = max(yy(:)) - min(yy(:));
    pad = max(0.12 * span, 1e-4);
    ylim(ax, [min(yy(:)) - pad, max(yy(:)) + pad]);
end

function styleAxes(ax, fontName)
    set(ax, 'FontName', fontName, 'FontSize', 11, 'LineWidth', 0.8);
    ax.XLabel.FontName = fontName;
    ax.YLabel.FontName = fontName;
    ax.Title.FontName = fontName;
    ax.Title.FontSize = 12;
    ax.Title.Interpreter = 'none';
    ax.XLabel.Interpreter = 'none';
    ax.YLabel.Interpreter = 'none';
    if ~isempty(ax.Legend)
        ax.Legend.FontName = fontName;
        ax.Legend.FontSize = 9;
        ax.Legend.Interpreter = 'none';
        ax.Legend.AutoUpdate = 'off';
        ax.Legend.Box = 'off';
    end
end

function saveFigure(fig, figDir, baseName)
    pngPath = fullfile(figDir, [baseName, '.png']);
    figPath = fullfile(figDir, [baseName, '.fig']);
    exportgraphics(fig, pngPath, 'Resolution', 200);
    savefig(fig, figPath);
end

function banner(fid, txt)
    bar = repmat('=', 1, 78);
    logf(fid, '\n%s\n%s\n%s\n', bar, txt, bar);
end

function printPair(fid, name, leftVal, rightVal, fmt)
    logf(fid, '%-22s  ', name);
    logf(fid, fmt, leftVal);
    logf(fid, '  ');
    logf(fid, fmt, rightVal);
    logf(fid, '\n');
end

function logf(fid, fmt, varargin)
%LOGF  同一段文字同時送到 Command Window 與 UTF-8 結果檔。
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end
