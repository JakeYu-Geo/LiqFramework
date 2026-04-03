%% ========== 主程序：ISA 模型模拟应力控制不排水循环扭剪（单位：kPa）==========
clear all; clc; close all;

% ---------- 模型参数（DEM，应力单位统一为 kPa）----------
props = [
    28*pi/180;  % 1  phic  临界摩擦角 (rad)
    0;          % 2  [已移除] 占位
    4000;       % 3  hs    颗粒硬度 (kPa)  原 2.6e6 Pa → 2600 kPa
    0.27;        % 4  n     指数
    0.61;       % 5  ed0   最小孔隙比 (p=0)
    0.89;       % 6  ec0   临界孔隙比 (p=0)
    1.00;       % 7  ei0   最大孔隙比 (p=0)
    0.2;       % 8  alpha 密度敏感指数
    4.0;        % 9  beta  密度敏感指数（用于fb）
    0.0;        % 10 mt    传统IS参数
    4.5;        % 11 mR    最大刚度放大因子
    1.0e-4;     % 12 R     屈服面直径（应变）
    0.08;        % 13 betaR0 硬化率
    0.9;        % 14 chi0  初始 chi
    3.7;       % 15 chi  chi
    0.02;       % 16 eaccPar 循环累积参数
    0;          % 17 
    300;        % 18 cz    Z 张量演化速率
    3;        % 19 beta_hor 各向同性硬化修正
    0.11;       % 20 gamma0_max 液化后最大剪应变
];
ntens = 6;      % Voigt 记号 6 分量

% ---------- 初始状态（有效应力，单位 kPa）----------
p0   = 100;                       % 初始有效围压 100 kPa
stress0 = [-p0; -p0; -p0; 0; 0; 0]; % 有效应力（拉正压负）
void0 = 0.693;                      % 初始孔隙比
stran = zeros(ntens,1);           % 总应变

% 状态变量向量（长度≥75）
statev = zeros(80,1);
statev(1) = void0;                % 孔隙比
statev(2) = 0;                   % 原pw位置，不再使用
statev(51:56) = 0;              % 初始总应变
statev(79) = 20;
statev(80) = 0;
statev = Initialasv(statev, props, ntens);  % 初始化 intergranular strain

% ---------- 应力控制循环参数 ----------
q_target = 20;                  % 目标偏应力幅值 (kPa)
axial_inc = 2e-4;              % 每步轴向应变增量（压缩为正）
total_steps = 15100;            % 总计算步数

% 初始化加载方向：+1 为压缩，-1 为拉伸
dir = 1;
hist = struct('p',[],'q',[],'e1',[],'s1',[],'cycle',[],'step',[]);

fprintf('开始应力控制循环加载（单位：kPa）...\n');
for step = 1:total_steps
    % 当前轴向应变增量
    deps1 = dir * axial_inc;
    % 不排水条件：deps_v = deps1 + 2*deps3 = 0
    deps3 = -0.5 * deps1;
    dstran = [0; 0; 0; -deps1; 0; 0];
    
    % 调用 IS 本构积分
    [stress, statev] = umatIS(stress0, statev, stran, dstran, props);
    
    % 更新状态
    stress0 = stress;
    stran = stran + dstran;
    
    % 计算当前偏应力 q
    [p, q] = pq(stress0);
    
    % --- 应力控制：检查是否达到目标值，反转方向 ---
    if dir == 1 && q >= q_target
        dir = -1;
        fprintf('  压缩半程结束，转向拉伸，步数 %d\n', step);
    elseif dir == -1 && q <= -q_target
        dir = 1;
        fprintf('  拉伸半程结束，转向压缩，步数 %d\n', step);
    end
    
    % 记录历史
    hist.p   = [hist.p;   p];
    hist.q   = [hist.q;   q];
    hist.e1  = [hist.e1;  -stran(4)];
    hist.s1  = [hist.s1;  stress0(4)];
    hist.step = [hist.step; step];
end

% ---------- 绘图 ----------
figure;
subplot(211); plot(-hist.e1 * 2, -hist.q); grid on;
xlabel('\epsilon_1'); ylabel('q (kPa)'); title('应力-应变');
subplot(212); plot(hist.p, hist.q); grid on; axis equal;
xlabel('p (kPa)'); ylabel('q (kPa)'); title('有效应力路径');

%% ========== 本地函数 ==========
% ------- ISA 本构积分核心-------
function [stress, statev] = umatIS(stress, statev, stran, dstran, props)
    % props 顺序见主程序注释，props(2) 已移除（占位0）
    % 所有应力单位：kPa；时间步长固定为 1（子步内 dt_sub = 1）
    persistent Mc c_ sq2 sq3 pcutmin Isym IsymES delta 
    if isempty(Mc)
        sq2 = sqrt(2); sq3 = sqrt(3); pcutmin = 0.1;   % kPa，张力截断
        phic = props(1); Mc = 6*sin(phic)/(3-sin(phic)); c_ = 3/(3+Mc);
        Isym = Isym1(); IsymES = Isym; IsymES(4:6,:) = 2*IsymES(4:6,:);
        delta = [1;1;1;0;0;0];
    end
    
    % 提取状态变量
    void = statev(1);
    hb = statev(3:8); cb = statev(9:14); zb = statev(15:20);
    eacc = statev(30); eps = stran;
    trS = trace(stress,3);devS = dev(stress);[p,q] = pq(stress);
    expB = exp(-(-trS/props(3))^props(4)); 
    ed = props(5)*expB; ec = props(6)*expB;
    rb = devS/p; rbiso = rb/(sq2/sq3*Mc);
    normrb = normS(rbiso);
    gS = getThetaS(devS, c_); gS = max(min(gS,1),c_);
    Fm = 1 + (normrb/gS)*(gS-1);
    fd0 = ((void-ed)/(ec-ed))^props(8);
    gamma0 = statev(77); epsv_mono = statev(76);Liqflag = statev(78);
    fe0 = (ec/void)^props(9);
    a = sq3*(3-sin(props(1)))/(2*sq2*sin(props(1)));
    trS = trace(stress,3);
    hatT = stress/trS;
    expB = exp(-(-trS/props(3))^props(4));
    ei = props(7)*expB;
    fb = (props(3)/props(4)*(props(7)/props(6))^props(9)*(1+ei)/ei * ...
                  (-trS/props(3))^(1-props(4))) / ...
                  (3 + a^2 - a*sq3*((props(7)-props(5))/(props(6)-props(5)))^props(8));
    % 张量 L 和 EE
    fs = fb*fe0 / normS(hatT)^2;
    Nb = unitE(hb-cb);
    rho = 1 - normE(props(12)*Nb - hb)/(2*props(12));
    rho = max(0, min(1, rho));
    pmin = statev(79);
    
    % 子步划分
    nDstran = normE(dstran);
    nsub = max(floor(nDstran/1e-5), 1);
    nsub = min(1, 5000);  % 不使用子步
    dstranU = dstran / nsub;
    DstranUU = unitE(dstranU);
    gamma_bar = stran - trace(stran) / 3 * delta;
    dt_sub = 1;          % 每个子步的时间固定为1（隐含每步 dt=1）
    
    for i = 1:nsub
        if -trace(stress,3)/3 < pcutmin 
            epsv_mono = epsv_mono + trace(dstranU);
            gamma_bar = gamma_bar + dstranU - trace(dstranU) / 3 * delta;
            normrb = normS(gamma_bar); gS = getThetaS(gamma_bar, c_);
            gS = max(min(gS,1),c_); g = 1 + (normrb/gS)*(gS-1);
            [~,tmp] = pq(gamma_bar); gamma_gen = abs(tmp) ;
            if ~Liqflag
                Liqflag = true;
                [~,tmp] = pq(gamma_bar);gamma_gen = abs(tmp) ;
                statev(80) = gamma_gen - gamma0 / g;
                gamma0 = gamma0 + 0.2 * mb(props(20) - gamma0);
                epsv_mono = 0;
            end            
            if sign(doubleSS((dstranU - trace(dstranU) / 3 * delta), (stran - trace(stran) / 3 * delta))) < 0
                gamma_gen = gamma_gen - statev(80);
            end
            if gamma_gen * g + epsv_mono * 1 > gamma0 
                stress = dev(stress) + (dstranU) * fs / 100 - delta * pcutmin;
                trS = trace(stress,3);
                hatT = stress/trS;
                EEhat = Fm^2*Isym + a^2 * (hatT*hatT');
                EE = EEhat * fs;
                Nhp = fd0*fs*a*Fm*(hatT + (hatT - delta/3));
                term2y = mb(doubleEE(Nb, DstranUU));
                chi_cur = props(14) + eacc*(props(15)-props(14));
                term3y = rho^chi_cur; yh = term2y * term3y;
                vN1 = Nhp; v1N = DstranUU'; v1N(4:6) = v1N(4:6)/2;
                Jmom = (props(11) + (1-props(11))*yh) * (EE + term3y * (vN1 * v1N));
                [dp,~] = pq(Jmom * dstranU);
                while dp < 0.05
                    stress = dev(stress) + (dstranU) * fs / 100 - delta * pcutmin;
                    trS = trace(stress,3);
                    hatT = stress/trS;
                    EEhat = Fm^2*Isym + a^2 * (hatT*hatT');
                    EE = EEhat * fs;
                    Nhp = fd0*fs*a*Fm*(hatT + (hatT - delta/3));
                    vN1 = Nhp;
                    Jmom = (props(11) + (1-props(11))*yh) * (EE + term3y * (vN1 * v1N));
                    [dp,~] = pq(Jmom * dstranU);
                end
            end
            
        else
            
            stress0 = stress; hb0 = hb; cb0 = cb;
            
            % --- 弹性刚度 & 特征状态 ---
            
            trS = trace(stress,3);devS = dev(stress);
            hatT = stress/trS;
            expB = exp(-(-trS/props(3))^props(4)); 
            ed = props(5)*expB; ec = props(6)*expB;
            gS = getThetaS(devS, c_); gS = max(min(gS,1),c_);
            [p,q] = pq(stress); 
            rb = devS/p; rbiso = rb/(sq2/sq3*Mc);
            normrb = normS(rbiso);
            % 材料常数 a, Fm
            a = sq3*(3-sin(props(1)))/(2*sq2*sin(props(1)));
            gS = getThetaS(devS, c_); gS = max(min(gS,1),c_);
            Fm = 1 + (normrb/gS)*(gS-1); Fm2 = Fm^2;
            % 特征孔隙比（Bauer 公式，应力单位 kPa）
             ei = props(7)*expB;
            % 压硬因子 fb（单位适应 kPa）
            fb = (props(3)/props(4)*(props(7)/props(6))^props(9)*(1+ei)/ei * ...
                  (-trS/props(3))^(1-props(4))) / ...
                  (3 + a^2 - a*sq3*((props(7)-props(5))/(props(6)-props(5)))^props(8));
            fe0 = (ec/void)^props(9); fd0 = ((void-ed)/(ec-ed))^props(8);
            % 循环移动修正
            Nb = unitE(hb-cb); termz0 = doubleEE(zb,Nb);
            termz = mb(-termz0); 
            fe = fe0 - mb(fe0-1)*termz; 
            fd = fd0 + mb(1-fd0)*termz;
            % 张量 L 和 EE
            fs = fb*fe / normS(hatT)^2;
            EEhat = Fm2*Isym + a^2 * (hatT*hatT');
            EE = EEhat * fs;
            Nhp = fd*fs*a*Fm*(hatT + (hatT - delta/3));
            
            % --- ISA 塑性（dt_sub=1）---
            hbtrial = hb + dstranU;
            if normE(hbtrial-cb) < props(12)/2   % 弹性
                hb = hb + dstranU;
                Load = true; yh = 0;
            else                                 % 塑性
                Load = false;
                Nb = unitE(hbtrial-cb); hbb = props(12)*Nb;
                vec = hbb - hb;
                if normE(vec)==0, vec = Nb; end
                betaR = props(19) + (props(13)-props(19)) * abs(doubleEE(unitE(vec),unitE(hb))) * (1-termz);
                % 更新 cb, hb（dt_sub = 1）
                cbmax = props(12)/2 * DstranUU;
                cbbar = (cbmax-cb)/props(12)*betaR;
                dotgamma = mb(doubleEE(Nb, dstranU)) / (1 + doubleEE(cbbar,Nb));  % dt_sub=1 隐含
                cb = cb0 + betaR/props(12)*(cbmax-cb0)*dotgamma*dt_sub / (1 + betaR/props(12)*dotgamma*dt_sub);
                hb = (hb0 + (dstranU - dotgamma*dt_sub/(props(12)/2)*(hb0-cb))) / (1 + dotgamma*dt_sub/(props(12)/2));
                if dotgamma>0
                    Nb = unitE(hb-cb); cb = hb - props(12)/2*Nb;
                end
                % 校正
                if dt_sub>0 && dotgamma>0 && (normE(cb)>props(12)/2-1e-6*props(12) || ...
                                              doubleEE(cb,DstranUU)>props(12)/2-1e-6*props(12) || ...
                                              normE(hb)/props(12)>1-1e-6)
                    cb = props(12)/2 * unitE(cb);
                    Nb = unitE(hb-cb); hb = cb + props(12)/2*Nb;
                end
                Nb = unitE(hb-cb);
            end
            
            % rho 因子
            if Load
                rho = 0;
            else
                rho = 1 - normE(props(12)*Nb - hb)/(2*props(12));
                rho = max(0, min(1, rho));
            end
            
            % --- 雅可比矩阵 ---
            if Load
                Jmom = props(11) * EE;
            else
                term2y = mb(doubleEE(Nb, DstranUU));
                chi_cur = props(14) + eacc*(props(15)-props(14));
                term3y = rho^chi_cur; yh = term2y * term3y;
                vN1 = Nhp; v1N = DstranUU'; v1N(4:6) = v1N(4:6)/2;
                Jmom = (props(11) + (1-props(11))*yh) * (EE + term3y * (vN1 * v1N));
            end
            
            % 应力更新
            stress = stress0 + Jmom * dstranU;
            ndU = normE(dstranU);
            eacc = eacc + (props(16)/props(12)) * ((1-yh)-eacc) * ndU;
            rat = q/p/(fd0*Mc*Fm) - 1;
            zb = zb + props(18)*mb(rat) * (1*Nb - zb) * ndU;
            normz = normE(zb);
            if normz > 1, zb = Nb; end

            [p,~] = pq(stress0);
            [dp,~] = pq(Jmom * dstranU);
            if p < pmin && dp > 0
                if pmin - p < 0.3 && pmin < 20
                    pcutmin = p;
                end
                pmin = p;
            end
            Liqflag = false;
        end
        
        % 更新孔隙比（不排水体积不变，但需根据体积应变更新）
        void = void + trace(dstranU,3)*(1+void);
        eps = eps + dstranU;
    end
    
    % 存储状态变量
    statev(1) = void; statev(2) = 0;  % pw 位置置零
    statev(3:8) = hb; statev(9:14) = cb; statev(15:20) = zb;
    statev(30) = eacc; statev(32) = normE(zb);
    statev(51:56) = eps;
    statev(76) = epsv_mono;
    statev(77) = gamma0;
    statev(78) = Liqflag;
    statev(79) = pmin;
    [p,q] = pq(stress); 
    statev(69)=p; statev(70)=q; statev(71)=q/p/(Mc*gS);
    statev(74)=void-ec; statev(75)=normE(hb)/props(12);
end

% ------- 初始化 intergranular strain -------
function statev = Initialasv(statev, props, ntens)
R = props(12);
hb = zeros(ntens,1); cb = zeros(ntens,1);
hb(1:3) = -0.57735*R; hb = hb / normE(hb) * R;
cb(1:3) = hb(1:3) + 0.57735*R/2;
statev(3:2+ntens) = hb; statev(9:8+ntens) = cb;
end

%% ---------- 张量运算辅助函数（Voigt 记号）----------
function tr = trace(A, ~), tr = A(1)+A(2)+A(3); end
function d = dev(A), d = A - trace(A,3)/3 * [1;1;1;0;0;0]; end
function n = normS(A), n = sqrt(sum(A(1:3).^2) + 2*sum(A(4:6).^2)); end
function n = normE(A), n = sqrt(sum(A(1:3).^2) + 0.5*sum(A(4:6).^2)); end
function u = unitS(A), n = normS(A); if n==0, u=zeros(size(A)); else u=A/n; end; end
function u = unitE(A), n = normE(A); if n==0, u=zeros(size(A)); else u=A/n; end; end
function r = doubleSS(A,B), r = sum(A(1:3).*B(1:3)) + 2*sum(A(4:6).*B(4:6)); end
function r = doubleEE(A,B), r = sum(A(1:3).*B(1:3)) + 0.5*sum(A(4:6).*B(4:6)); end
function M = dyadSS(A,B,~), M = A*B'; end
function [p,q] = pq(T), d = dev(T); p = -trace(T,3)/3; q = -T(4); end
function y = mb(x), y = max(0,x); end

function Isym = Isym1()
Isym = diag([1,1,1,0.5,0.5,0.5]);
Isym(1:3,1:3) = 1/3;
Isym(1:3,1:3) = Isym(1:3,1:3) + eye(3)/2;
end

function g = getThetaS(A, c)
sq6 = sqrt(6);
B = [A(1:3); A(4:6)/2]; B = dev(B); B = unitE(B);
cos3t = B(1)^3 + B(2)^3 + B(3)^3 ...
      + 3*B(3)*B(5)^2 + 3*B(1)*(B(4)^2+B(5)^2) ...
      + 6*B(4)*B(5)*B(6) + 3*B(3)*B(6)^2 + 3*B(2)*(B(4)^2+B(6)^2);
cos3t = -sq6 * cos3t;
g = 2*c / ((1+c) - (1-c)*cos3t);
g = max(min(g,1), c);
end
