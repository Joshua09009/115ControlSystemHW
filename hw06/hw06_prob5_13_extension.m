%% ========================================================================
%  hw06_prob5_13_extension.m
%  Control Systems HW06 - Problem 5.13 (Unit 5A/5B/5C: Root Locus)
%  原題標準解答重製 + 進階工程延伸 (Baseline Verification + Extensions)
%
%  System (unity negative feedback):
%     D(s) = K (s+1)/(s+13)                 Lead compensator
%     G(s) = (s^2+81) / [s^2 (s^2+100)]     Flexible plant
%     L(s) = D(s)G(s),  T = L/(1+L),  M = U/R = D/(1+DG)
%
%  Modules
%     M1  原題精確求解與模態拆解 (Parts a-d)          -> Figure 1
%     M2  極零點翻轉強健性分析 (Pole-Zero Flip)       -> Figure 2
%     M3  控制力與 Notch Pre-filter 改良設計           -> Figure 3
%     M4  Markdown 格式總表 (Formatted Summary Table)
%
%  需求: Control System Toolbox, MATLAB R2016b+ (script 內 local function)
%  備註: 圖上文字與 Command Window 輸出使用英文，避免中文字型/編碼亂碼。
%% ========================================================================
clear; close all; clc;

%% ---------------- 0. 全域設定 (Global Settings) -------------------------
FS   = 12;          % 字體大小 FontSize
LW   = 1.5;         % 線寬 LineWidth
K_tb = 35;          % 課本解答增益 (textbook gain)
Kscan = linspace(0.1, 200, 50000);   % 高密度增益掃描 (high-density K sweep)
SAVE_FIG = true;    % 是否輸出 PNG 檔

set(groot, 'defaultAxesFontSize', FS, 'defaultLineLineWidth', LW, ...
           'defaultAxesLineWidth', 1.0, 'defaultTextFontSize', FS);

% ---- 名目系統多項式 (Nominal polynomials) ----
numD0 = [1 1];                 denD = [1 13];          % D(s)/K
numG  = [1 0 81];              denG = conv([1 0 0], [1 0 100]);
numL  = conv(numD0, numG);     denL = conv(denD, denG); % L(s)/K
D0 = tf(numD0, denD);   G = tf(numG, denG);

fprintf('=====================================================================\n');
fprintf(' HW06 Problem 5.13  --  Root Locus Baseline + Extensions\n');
fprintf(' MATLAB %s\n', version);
fprintf('=====================================================================\n\n');

%% ========================================================================
%  MODULE 1 : Baseline Verification (Parts a, b, c, d)
%% ========================================================================
fprintf('#####################  MODULE 1 : BASELINE  #####################\n\n');

z = roots(numL);  p = roots(denL);
n = numel(p);     m = numel(z);

% ---- (a) 漸近線 Asymptotes ----
sigma_a = (sum(p) - sum(z)) / (n - m);
phi_a   = (180 + 360*(0:(n-m-1))) / (n - m);
fprintf('[M1-a] Open-loop poles : %s\n', vec2str(p));
fprintf('[M1-a] Open-loop zeros : %s\n', vec2str(z));
fprintf('[M1-a] n - m = %d\n', n - m);
fprintf('[M1-a] Asymptote centroid sigma_a = %.6f\n', real(sigma_a));
fprintf('[M1-a] Asymptote angles phi_i     = %s deg\n', mat2str(phi_a, 6));

% ---- (a) 出發角 Angle of departure at s = +j10 ----
idx_j10 = find(abs(p - 10j) < 1e-6, 1);
[phi_dep, ang_detail] = departure_angle(z, p, idx_j10);
phi_dep_num = numeric_departure(numL, denL, p(idx_j10));
fprintf('[M1-a] Departure angle at +j10 (angle condition) = %.4f deg\n', phi_dep);
fprintf('       sum(angle to zeros) = %.4f deg, sum(angle to other poles) = %.4f deg\n', ...
        ang_detail(1), ang_detail(2));
fprintf('[M1-a] Departure angle at +j10 (numerical, K=1e-7) = %.4f deg\n', phi_dep_num);

% 到達角 Arrival angle at +j9 (補充)
idx_z9  = find(abs(z - 9j) < 1e-6, 1);
psi_arr = arrival_angle(z, p, idx_z9);
fprintf('[M1-a] Arrival angle at zero +j9 (angle condition) = %.4f deg\n\n', psi_arr);

% ---- (b) 高密度掃描與分支追蹤 Branch tracking ----
R    = track_roots(numL, denL, Kscan);         % 5 x N, 每列為一條連續分支
zeta = -real(R) ./ abs(R);
R0   = R(:,1);
iRes = find(abs(abs(R0) - 10) < 0.5);          % 由 +-j10 出發 -> 諧振極點組
iDom = find(abs(R0) < 0.5);                    % 由原點出發   -> 主導極點組
iRig = setdiff(1:n, [iRes; iDom]);             % 由 -13 出發

zeta_res = min(zeta(iRes,:), [], 1);
zeta_dom = min(zeta(iDom,:), [], 1);
zeta_sys = min(zeta, [], 1);                   % 全部 5 根中最差的阻尼

[zr_max, ir] = max(zeta_res);
[zd_max, ~]  = max(zeta_dom);
id_first     = find(zeta_dom >= zd_max - 1e-9, 1);
[zs_max, is] = max(zeta_sys);

fprintf('[M1-b] Branch classification at K=%.1f: dominant=%s, resonant=%s, from -13=%s\n', ...
        Kscan(1), mat2str(iDom'), mat2str(iRes'), mat2str(iRig));
fprintf('[M1-b] Dominant pair : zeta_max = %.4f, first reached at K = %.4f\n', ...
        zd_max, Kscan(id_first));
kk = Kscan(zeta_dom > 0.5);
if ~isempty(kk)
    fprintf('[M1-b] Dominant pair : zeta > 0.5 for K in [%.4f, %.4f] (within sweep)\n', kk(1), kk(end));
end
fprintf('[M1-b] Resonant pair : zeta_max = %.4f at K = %.4f, pole = %s\n', ...
        zr_max, Kscan(ir), cplx2str(R(iRes(1),ir)));
fprintf('[M1-b] ALL 5 roots   : max_K( min zeta ) = %.4f at K = %.4f\n', zs_max, Kscan(is));

% 延伸掃描 (extended sweep) 驗證 K -> infinity 也達不到 0.5
Kext   = logspace(-2, 5, 40000);
Rext   = track_roots(numL, denL, Kext);
iResE  = find(abs(abs(Rext(:,1)) - 10) < 0.5);
zeExt  = min(-real(Rext(iResE,:)) ./ abs(Rext(iResE,:)), [], 1);
[zeE_max, ieE] = max(zeExt);
fprintf('[M1-b] Extended sweep K in [1e-2, 1e5]: resonant zeta_max = %.4f at K = %.3f\n', ...
        zeE_max, Kext(ieE));
fprintf('[M1-b] Resonant pole at K = 1e5: %s  (-> zero at +j9)\n\n', ...
        cplx2str(Rext(iResE(imag(Rext(iResE,end))>0), end)));

% ---- (c) 精確 K_opt (zeta_dom = sqrt(2)/2) ----
zt = sqrt(2)/2;
i0 = find(zeta_dom >= zt, 1);
K_lin = interp1(zeta_dom(i0-1:i0), Kscan(i0-1:i0), zt);
K_opt = fzero(@(K) zeta_dominant(K, numL, denL) - zt, [Kscan(i0-1), Kscan(i0)]);
fprintf('[M1-c] K for zeta_dom = 0.7071 : sweep-interp = %.5f, fzero = %.6f\n', K_lin, K_opt);
fprintf('[M1-c] zeta_dom at textbook K = 35 : %.5f\n\n', zeta_dominant(K_tb, numL, denL));

% ---- (d) 兩組增益的完整閉迴路分析 ----
res_K35  = analyze_cl(K_tb,  numD0, denD, numG, denG, sprintf('Nominal, K = %g (textbook)', K_tb));
res_Kopt = analyze_cl(K_opt, numD0, denD, numG, denG, sprintf('Nominal, K_opt = %.4f', K_opt));

% ---- Figure 1 : 全局 + 局部放大根軌跡 ----
Kplot = unique([logspace(-5,-1,600), Kscan, logspace(log10(200),5,4000)]);
Rp    = track_roots(numL, denL, Kplot);
pK35  = res_K35.poles;   pKopt = res_Kopt.poles;

fig1 = figure('Name','Fig1 Root Locus','Color','w','Position',[60 80 1300 560]);
cols = lines(n);
for sp = 1:2
    ax = subplot(1,2,sp); hold(ax,'on'); box(ax,'on');
    hb = gobjects(n,1);
    for b = 1:n
        hb(b) = plot(ax, real(Rp(b,:)), imag(Rp(b,:)), '-', 'Color', cols(b,:), 'LineWidth', LW);
    end
    hp  = plot(ax, real(p), imag(p), 'kx', 'MarkerSize', 11, 'LineWidth', 2);
    hz  = plot(ax, real(z), imag(z), 'ko', 'MarkerSize', 9,  'LineWidth', 1.5);
    hk1 = plot(ax, real(pKopt), imag(pKopt), 'rs', 'MarkerSize', 9, 'MarkerFaceColor','r');
    hk2 = plot(ax, real(pK35),  imag(pK35),  'bd', 'MarkerSize', 8, 'LineWidth', 1.5);
    if sp == 1
        sgrid([0.5 0.707], []);
        % 漸近線 asymptotes
        ha = plot(ax, real(sigma_a)*[1 1], [-30 30], 'k--', 'LineWidth', 1.0);
        axis(ax, [-20 4 -16 16]);
        title(ax, {'(a) Global root locus', 'L(s) = K(s+1)(s^2+81) / [s^2(s+13)(s^2+100)]'}, 'FontSize', FS-1);
        legend(ax, [hp hz hk1 hk2 ha], {'OL poles','OL zeros', ...
               sprintf('CL poles, K_{opt}=%.2f', K_opt), sprintf('CL poles, K=%g', K_tb), ...
               sprintf('Asymptotes, \\sigma_a=%.1f', real(sigma_a))}, 'Location','southwest', 'AutoUpdate','off');
    else
        % 諧振區局部放大 Resonance zoom
        plot(ax, real(R(iRes(imag(R(iRes,ir))>0), ir)), imag(R(iRes(imag(R(iRes,ir))>0), ir)), ...
             'mp', 'MarkerSize', 13, 'MarkerFaceColor', 'm');
        th = linspace(pi/2, pi, 200);
        for zz = [0.1 0.2 0.5]   % 參考阻尼射線 reference zeta rays
            ang = pi - acos(zz);
            plot(ax, [0 15*cos(ang)], [0 15*sin(ang)], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            text(ax, 11.2*cos(ang), 11.2*sin(ang), sprintf('\\zeta=%.1f', zz), 'Color', [0.4 0.4 0.4], 'FontSize', FS-2);
        end
        axis(ax, [-3 1 7.5 11.5]);
        title(ax, '(b) Zoom: resonance branch +j10 \rightarrow +j9', 'FontSize', FS-1);
        text(ax, -2.9, 7.75, sprintf('\\phi_{dep}(+j10) = %.2f^\\circ,  \\zeta_{res,max} = %.3f @ K = %.1f', ...
             phi_dep, zr_max, Kscan(ir)), 'FontSize', FS-2);
    end
    xline(ax, 0, 'k-', 'LineWidth', 0.8);
    xlabel(ax, 'Real axis (s^{-1})'); ylabel(ax, 'Imaginary axis (s^{-1})'); grid(ax,'on');
end
if SAVE_FIG, save_fig(fig1, 'Fig1_RootLocus'); end

%% ========================================================================
%  MODULE 2 : Pole-Zero Flip Robustness (Extension 1)
%% ========================================================================
fprintf('\n#####################  MODULE 2 : POLE-ZERO FLIP  #####################\n\n');
wz2_list  = [81 96 105];
case_name = {'Case 1: Nominal collocated (wz^2=81)', ...
             'Case 2: Near-cancellation (wz^2=96)', ...
             'Case 3: Pole-zero flip (wz^2=105)'};
T_end2 = 60;                       % 長時間模擬以顯示 Case 2 慢衰減 / Case 3 發散
t2 = 0:2e-3:T_end2;
res_case = cell(1,3);
Kz = unique([logspace(-5,-1,400), linspace(0.1,200,20000), logspace(log10(200),4,1500)]);

fig2 = figure('Name','Fig2 Pole-Zero Flip','Color','w','Position',[80 60 1250 900]);
cc = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188];
ystep = zeros(numel(t2),3);
for c = 1:3
    wz2   = wz2_list(c);
    numGc = [1 0 wz2];
    numLc = conv(numD0, numGc);
    zc = roots(numLc);  pc = roots(denL);
    iP = find(abs(pc - 10j) < 1e-6, 1);
    phi_c = departure_angle(zc, pc, iP);
    phi_c_num = numeric_departure(numLc, denL, pc(iP));
    fprintf('[M2] %s : wz = %.4f rad/s, wp = 10 rad/s\n', case_name{c}, sqrt(wz2));
    fprintf('     Departure angle at +j10 : analytic = %.4f deg, numerical = %.4f deg\n', phi_c, phi_c_num);

    rc = analyze_cl(K_tb, numD0, denD, numGc, denG, sprintf('%s, K=%g', case_name{c}, K_tb));
    rc.wz2 = wz2;  rc.phi_dep = phi_c;
    res_case{c} = rc;

    % 局部根軌跡 Zoomed locus
    Rc = track_roots(numLc, denL, Kz);
    ax = subplot(2,2,c); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for b = 1:n
        plot(ax, real(Rc(b,:)), imag(Rc(b,:)), '-', 'Color', cc(c,:), 'LineWidth', LW);
    end
    plot(ax, real(pc), imag(pc), 'kx', 'MarkerSize', 11, 'LineWidth', 2);
    plot(ax, real(zc), imag(zc), 'ko', 'MarkerSize', 9,  'LineWidth', 1.5);
    plot(ax, real(rc.poles), imag(rc.poles), 'rs', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
    % 出發方向箭頭 departure direction
    quiver(ax, 0, 10, 0.6*cosd(phi_c), 0.6*sind(phi_c), 0, 'Color', 'm', 'LineWidth', 2, 'MaxHeadSize', 0.8);
    xline(ax, 0, 'k-', 'LineWidth', 1.0);
    axis(ax, [-2 0.6 8.5 10.8]);
    xlabel(ax, 'Real axis (s^{-1})'); ylabel(ax, 'Imaginary axis (s^{-1})');
    title(ax, sprintf('%s\n\\phi_{dep} = %.1f^\\circ, Re\\{s_{res}\\}(K=%g) = %+.4f', ...
          case_name{c}, phi_c, K_tb, real(rc.p_res)), 'FontSize', FS-1);

    Tc = feedback(K_tb*D0*tf(numGc, denG), 1);
    ystep(:,c) = step(Tc, t2);
end
ax = subplot(2,2,4); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
hs = gobjects(3,1);
for c = 1:3
    hs(c) = plot(ax, t2, ystep(:,c), '-', 'Color', cc(c,:), 'LineWidth', LW);
end
yline(ax, 1, 'k:', 'LineWidth', 1);
xlabel(ax, 'Time (s)'); ylabel(ax, 'y(t)');
title(ax, sprintf('Step response at K = %g (0 - %g s)', K_tb, T_end2), 'FontSize', FS-1);
legend(ax, hs, {'Case 1 (\omega_z=9)','Case 2 (\omega_z=9.80)','Case 3 (\omega_z=10.25)'}, 'Location','northwest', 'AutoUpdate','off');
ylim(ax, [min(-0.2, min(ystep(:))-0.1), max(1.6, max(ystep(:))+0.1)]);
if SAVE_FIG, save_fig(fig2, 'Fig2_PoleZeroFlip'); end

% 發散/衰減包絡量化 Envelope quantification
fprintf('\n[M2] Residual envelope |y-1| in windows (max):\n');
fprintf('     %-38s  [5,10]s     [25,30]s    [55,60]s\n', 'Case');
for c = 1:3
    w1 = t2>=5 & t2<=10; w2 = t2>=25 & t2<=30; w3 = t2>=55 & t2<=60;
    e = abs(ystep(:,c)-1);
    fprintf('     %-38s  %.4e  %.4e  %.4e\n', case_name{c}, max(e(w1)), max(e(w2)), max(e(w3)));
end

%% ========================================================================
%  MODULE 3 : Control Effort & Notch Pre-filter Redesign (Extension 2)
%% ========================================================================
fprintf('\n#####################  MODULE 3 : NOTCH PRE-FILTER  #####################\n\n');
t3 = 0:1e-3:6;

% Design A : Baseline
T_A = feedback(K_tb*D0*G, 1);
M_A = feedback(K_tb*D0, G);             % U/R = D/(1+DG)

% 對準閉迴路諧振極點 Target the closed-loop resonant pole of Design A
pA   = pole(T_A);
p_r  = pick_resonant(pA);
w0   = abs(p_r);
zr   = -real(p_r)/w0;
zeta_p = 0.7;                            % 濾波器分母阻尼 (well-damped filter poles)
F = tf([1 2*zr*w0 w0^2], [1 2*zeta_p*w0 w0^2]);
fprintf('[M3] Closed-loop resonant pole (Design A) p_r = %s\n', cplx2str(p_r));
fprintf('[M3] Pre-filter F(s) = (s^2 + 2*%.5f*%.5f s + %.5f^2) / (s^2 + 2*%.2f*%.5f s + %.5f^2)\n', ...
        zr, w0, w0, zeta_p, w0, w0);
fprintf('     F(0) = %.4f (unity DC gain), F(inf) = 1\n', dcgain(F));
fprintf('     Pre-filter is outside the loop -> closed-loop poles / stability unchanged.\n');

% Design B : Lead + Notch pre-filter (nominal plant)
T_B = T_A * F;     M_B = M_A * F;
% Design B under plant drift (wz^2 = 96): 預濾波器未重調 (filter not re-tuned)
Gd   = tf([1 0 96], denG);
T_Bd = feedback(K_tb*D0*Gd, 1) * F;   M_Bd = feedback(K_tb*D0, Gd) * F;

yA = step(T_A, t3);  uA = step(M_A, t3);
yB = step(T_B, t3);  uB = step(M_B, t3);
yBd = step(T_Bd, t3); uBd = step(M_Bd, t3);

mA  = effort_metrics('Design A (Baseline)',            T_A,  M_A,  t3, yA,  uA,  p_r);
mB  = effort_metrics('Design B (Lead + Notch F(s))',   T_B,  M_B,  t3, yB,  uB,  p_r);
mBd = effort_metrics('Design B on drifted plant (96)', T_Bd, M_Bd, t3, yBd, uBd, pick_resonant(pole(feedback(K_tb*D0*Gd,1))));

fig3 = figure('Name','Fig3 Notch Prefilter','Color','w','Position',[100 60 1000 820]);
ax = subplot(2,1,1); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, t3, yA, '-',  'Color', [0 0.447 0.741], 'LineWidth', LW);
plot(ax, t3, yB, '-',  'Color', [0.85 0.325 0.098], 'LineWidth', LW);
plot(ax, t3, yBd,'--', 'Color', [0.466 0.674 0.188], 'LineWidth', LW);
yline(ax, 1, 'k:', 'LineWidth', 1);
xlabel(ax, 'Time (s)'); ylabel(ax, 'y(t)');
title(ax, sprintf('Step output: Baseline vs Notch pre-filter (\\omega_0 = %.3f rad/s)', w0), 'FontSize', FS);
legend(ax, {sprintf('A: Baseline (A_{res} = %.4f)', mA.A_res), ...
            sprintf('B: Lead + Notch F(s) (A_{res} = %.2e)', mB.A_res), ...
            sprintf('B on drifted plant \\omega_z^2=96 (A_{res} = %.4f)', mBd.A_res)}, 'Location','southeast');

ax = subplot(2,1,2); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, t3, uA, '-',  'Color', [0 0.447 0.741], 'LineWidth', LW);
plot(ax, t3, uB, '-',  'Color', [0.85 0.325 0.098], 'LineWidth', LW);
plot(ax, t3, uBd,'--', 'Color', [0.466 0.674 0.188], 'LineWidth', LW);
xlabel(ax, 'Time (s)'); ylabel(ax, 'u(t)');
title(ax, 'Control effort u(t) = M(s)F(s)R(s)', 'FontSize', FS);
legend(ax, {sprintf('A: u(0^+)=%.2f, RMS_{[0,6]}=%.3f, RMS_{[3,6]}=%.4f', mA.u0, mA.uRMS, mA.uRMS_tail), ...
            sprintf('B: u(0^+)=%.2f, RMS_{[0,6]}=%.3f, RMS_{[3,6]}=%.4f', mB.u0, mB.uRMS, mB.uRMS_tail), ...
            sprintf('B drift: u(0^+)=%.2f, RMS_{[0,6]}=%.3f, RMS_{[3,6]}=%.4f', mBd.u0, mBd.uRMS, mBd.uRMS_tail)}, ...
       'Location','northeast');
% 插入初始段放大 inset of initial spike
axi = axes('Position',[0.29 0.25 0.17 0.11]); hold(axi,'on'); box(axi,'on'); grid(axi,'on');
plot(axi, t3, uA, '-', 'Color', [0 0.447 0.741]); plot(axi, t3, uB, '-', 'Color', [0.85 0.325 0.098]);
xlim(axi, [0 0.5]); title(axi, 'zoom 0-0.5 s', 'FontSize', FS-3); set(axi, 'FontSize', FS-3);
if SAVE_FIG, save_fig(fig3, 'Fig3_NotchPrefilter'); end

%% ========================================================================
%  MODULE 4 : Formatted Summary Table (Markdown)
%% ========================================================================
fprintf('\n#####################  MODULE 4 : SUMMARY TABLES (copy below)  #####################\n\n');
fprintf('<!-- BEGIN SUMMARY -->\n');
fprintf('**Table 1. Key root-locus quantities (nominal)**\n\n');
fprintf('| Quantity | Value |\n|---|---|\n');
fprintf('| sigma_a | %.4f |\n', real(sigma_a));
fprintf('| phi_i (deg) | %s |\n', mat2str(phi_a, 5));
fprintf('| phi_dep at +j10 (deg) | %.4f |\n', phi_dep);
fprintf('| psi_arr at +j9 (deg) | %.4f |\n', psi_arr);
fprintf('| Dominant zeta_max / first K | %.4f / %.4f |\n', zd_max, Kscan(id_first));
fprintf('| Resonant zeta_max / K | %.4f / %.4f |\n', zr_max, Kscan(ir));
fprintf('| Resonant zeta_max, K in [1e-2,1e5] | %.4f |\n', zeE_max);
fprintf('| max_K min(zeta over 5 roots) / K | %.4f / %.4f |\n', zs_max, Kscan(is));
fprintf('| K_opt (zeta_dom = 0.7071) | %.6f |\n', K_opt);
fprintf('| zeta_dom at K = 35 | %.5f |\n\n', zeta_dominant(K_tb, numL, denL));

fprintf('**Table 2. Closed-loop poles and step metrics**\n\n');
fprintf('| Experiment | K | CL poles | zeta_dom | wn_dom | zeta_res | wn_res | %%OS | t_r (s) | t_s 2%% (s) | A_res (10 rad/s mode) |\n');
fprintf('|---|---|---|---|---|---|---|---|---|---|---|\n');
all_rows = [{res_K35, res_Kopt}, res_case];
for k = 1:numel(all_rows)
    r = all_rows{k};
    fprintf('| %s | %.4f | %s | %.4f | %.4f | %.4f | %.4f | %s | %s | %s | %.4e |\n', ...
        r.label, r.K, vec2str(r.poles), r.zeta_dom, r.wn_dom, r.zeta_res, r.wn_res, ...
        num_or_na(r.OS,'%.2f'), num_or_na(r.tr,'%.4f'), num_or_na(r.ts,'%.4f'), r.A_res);
end

fprintf('\n**Table 3. Pole-zero flip (K = %g)**\n\n', K_tb);
fprintf('| Case | wz | phi_dep (deg) | Re{s_res} | Im{s_res} | Stable? |\n|---|---|---|---|---|---|\n');
for c = 1:3
    r = res_case{c};
    fprintf('| %s | %.4f | %.4f | %+.5f | %.5f | %s |\n', case_name{c}, sqrt(r.wz2), r.phi_dep, ...
            real(r.p_res), imag(r.p_res), tern(r.stable,'Yes','NO'));
end

fprintf('\n**Table 4. Control effort & residual vibration (0-6 s)**\n\n');
fprintf('| Design | %%OS | t_s 2%% (s) | A_res (modal) | max|y-1| on [4,6]s | u(0+) | max|u| | u RMS [0,6] | u RMS [3,6] |\n');
fprintf('|---|---|---|---|---|---|---|---|---|\n');
for mm = {mA, mB, mBd}
    q = mm{1};
    fprintf('| %s | %s | %s | %.4e | %.4e | %.4f | %.4f | %.4f | %.4e |\n', q.label, ...
        num_or_na(q.OS,'%.2f'), num_or_na(q.ts,'%.4f'), q.A_res, q.res_tail, q.u0, q.umax, q.uRMS, q.uRMS_tail);
end
fprintf('\n<!-- END SUMMARY -->\n');
fprintf('\nDone. Figures saved: %d\n', SAVE_FIG);

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================
function R = track_roots(num, den, Kvec)
% 依連續性追蹤根軌跡分支 (continuity-based branch tracking)
    n = numel(den) - 1;
    numP = [zeros(1, numel(den)-numel(num)), num];
    P = perms(1:n);
    R = zeros(n, numel(Kvec));
    prev = roots(den + Kvec(1)*numP);
    R(:,1) = prev;
    for k = 2:numel(Kvec)
        r = roots(den + Kvec(k)*numP);
        cost = sum(abs(r(P) - repmat(prev.', size(P,1), 1)), 2);
        [~, ib] = min(cost);
        r = r(P(ib,:));
        R(:,k) = r;  prev = r;
    end
end

function [phi, detail] = departure_angle(z, p, idx)
% 角度條件求出發角 (angle condition): phi = 180 + sum(angles to zeros) - sum(angles to other poles)
    p0 = p(idx);  po = p;  po(idx) = [];
    sz = sum(rad2deg(angle(p0 - z)));
    sp = sum(rad2deg(angle(p0 - po)));
    phi = wrap180(180 + sz - sp);
    detail = [sz, sp];
end

function psi = arrival_angle(z, p, idx)
% 到達角 (angle of arrival) at zero z(idx)
    z0 = z(idx);  zo = z;  zo(idx) = [];
    psi = wrap180(180 + sum(rad2deg(angle(z0 - p))) - sum(rad2deg(angle(z0 - zo))));
end

function phi = numeric_departure(num, den, p0)
% 以極小 K 數值驗證出發方向 (numerical check with tiny K)
    numP = [zeros(1, numel(den)-numel(num)), num];
    Ks = 1e-7;
    r = roots(den + Ks*numP);
    [~, k] = min(abs(r - p0));
    phi = rad2deg(angle(r(k) - p0));
end

function a = wrap180(a)
    a = mod(a + 180, 360) - 180;
end

function zd = zeta_dominant(K, numL, denL)
% 主導極點阻尼比: |s|<5 的複數極點 (valid for K ~ 1-50)
    numP = [zeros(1, numel(denL)-numel(numL)), numL];
    r = roots(denL + K*numP);
    rc = r(abs(r) < 5 & imag(r) > 1e-9);
    if isempty(rc), zd = 1; else, zd = min(-real(rc)./abs(rc)); end
end

function pr = pick_resonant(pl)
% 諧振極點: 上半平面中 |Im| > 5 且最接近虛軸者
    c = pl(imag(pl) > 5);
    [~, k] = max(real(c));
    pr = c(k);
end

function out = analyze_cl(K, numD0, denD, numG, denG, label)
% 閉迴路極點、模態、步階性能與殘餘震盪 (CL poles, modes, step metrics)
    D = K*tf(numD0, denD);   G = tf(numG, denG);
    T = feedback(D*G, 1);
    pl = pole(T);
    [~, ord] = sort(abs(pl));  pl = pl(ord);
    out.label = label;  out.K = K;  out.poles = pl;
    out.stable = all(real(pl) < 0);

    fprintf('---- %s ----\n', label);
    fprintf('   %-24s %10s %10s   %s\n', 'Pole', 'zeta', 'wn(rad/s)', 'Mode');
    for k = 1:numel(pl)
        zk = -real(pl(k))/abs(pl(k));
        if imag(pl(k)) > 5 || imag(pl(k)) < -5, md = 'resonant';
        elseif abs(imag(pl(k))) > 1e-9, md = 'dominant';
        else, md = 'real'; end
        fprintf('   %-24s %10.5f %10.5f   %s\n', cplx2str(pl(k)), zk, abs(pl(k)), md);
    end
    pr = pick_resonant(pl);
    out.p_res = pr;  out.zeta_res = -real(pr)/abs(pr);  out.wn_res = abs(pr);
    pd = pl(abs(pl) < 5 & imag(pl) > 1e-9);
    if isempty(pd)   % 主導極點為實根 (real dominant)
        pdr = pl(abs(imag(pl)) < 1e-9);  [~, k] = max(real(pdr)); pd = pdr(k);
    end
    out.zeta_dom = -real(pd(1))/abs(pd(1));  out.wn_dom = abs(pd(1));

    if out.stable
        S = stepinfo(T, 'SettlingTimeThreshold', 0.02, 'RiseTimeLimits', [0.1 0.9]);
        out.OS = S.Overshoot;  out.tr = S.RiseTime;  out.ts = S.SettlingTime;
        out.Tpeak = S.PeakTime; out.Ypeak = S.Peak;
    else
        out.OS = NaN; out.tr = NaN; out.ts = NaN; out.Tpeak = NaN; out.Ypeak = NaN;
    end
    out.A_res = modal_amp(T, pr);
    fprintf('   Stable: %s | %%OS = %s | t_r = %s s | t_s(2%%) = %s s | t_peak = %s s\n', ...
        tern(out.stable,'Yes','NO'), num_or_na(out.OS,'%.3f'), num_or_na(out.tr,'%.4f'), ...
        num_or_na(out.ts,'%.4f'), num_or_na(out.Tpeak,'%.4f'));
    fprintf('   Resonant-mode amplitude in y(t) (2|residue|) A_res = %.5e, freq = %.4f rad/s\n\n', ...
        out.A_res, imag(pr));
end

function A = modal_amp(sys, p_target)
% 步階響應中某複數模態的振幅 2|Res{Y(s), p}|
    [num, den] = tfdata(sys * tf(1, [1 0]), 'v');
    pl = roots(den);
    [~, k] = min(abs(pl - p_target));
    pk = pl(k);
    r = polyval(num, pk) / polyval(polyder(den), pk);
    A = 2*abs(r);
end

function m = effort_metrics(label, T, M, t, y, u, p_r)
    m.label = label;
    S = stepinfo(y, t, 1, 'SettlingTimeThreshold', 0.02);
    m.OS = S.Overshoot;  m.ts = S.SettlingTime;
    m.A_res = modal_amp(T, p_r);
    m.res_tail = max(abs(y(t >= 4) - 1));
    m.u0 = u(1);   m.umax = max(abs(u));
    m.uRMS = sqrt(mean(u.^2));
    m.uRMS_tail = sqrt(mean(u(t >= 3).^2));
    m.u_ss = dcgain(M);
    fprintf('[M3] %-32s %%OS=%7.3f  t_s=%7.4f s  A_res=%.4e  max|y-1|[4,6]=%.4e\n', ...
            label, m.OS, m.ts, m.A_res, m.res_tail);
    fprintf('     %-32s u(0+)=%8.4f  max|u|=%8.4f  RMS[0,6]=%8.4f  RMS[3,6]=%.4e  u_ss=%.3g\n', ...
            '', m.u0, m.umax, m.uRMS, m.uRMS_tail, m.u_ss);
end

function s = cplx2str(c)
    if abs(imag(c)) < 1e-9, s = sprintf('%.4f', real(c));
    else, s = sprintf('%.4f %s %.4fj', real(c), char(43 + 2*(imag(c)<0)), abs(imag(c))); end
end

function s = vec2str(v)
% 共軛對合併顯示 (compact, conjugate pairs as a +- bj)
    v = v(:);  out = {};  used = false(size(v));
    for k = 1:numel(v)
        if used(k), continue; end
        if abs(imag(v(k))) < 1e-9
            out{end+1} = sprintf('%.4f', real(v(k))); used(k) = true; %#ok<AGROW>
        else
            out{end+1} = sprintf('%.4f+/-%.4fj', real(v(k)), abs(imag(v(k)))); used(k) = true; %#ok<AGROW>
            [~, j] = min(abs(v - conj(v(k))) + 1e9*used);  used(j) = true;
        end
    end
    s = strjoin(out, ', ');
end

function s = num_or_na(x, fmt)
    if isnan(x) || isinf(x), s = 'N/A'; else, s = sprintf(fmt, x); end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

function save_fig(fig, name)
    fdir = fullfile(pwd, 'figures');
    if ~isfolder(fdir), mkdir(fdir); end
    name = fullfile(fdir, name);
    try
        exportgraphics(fig, [name '.png'], 'Resolution', 300);
    catch
        print(fig, name, '-dpng', '-r300');
    end
end
