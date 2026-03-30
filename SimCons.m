clc; clear;

G0 = 125.0;
nu = 0.05;
M = 1.25;
c = 0.712;
lambdaC = 0.019;
e0 = 0.934;
xi = 0.7;
m = 0.01;
h0 = 7.05;
ch = 0.968;
nb = 1.1;
A0 = 0.704;
nd = 3.5;
zmax = 4;
% eA = -0.18;
cz = 600;
% cF = 5;
pAtmos = 101.4;
pMin = 2;
% qMax = 100;

% 初始化变量
s = {zeros(3, 3)};
p = [100.0]';
z = {zeros(3, 3)};
% F = {diag([2, -1, -1]) / sqrt(6) * 0.4};
e = {zeros(3, 3)};
sigma = {zeros(3, 3)};
alpha = zeros(3, 3);
alphaInit = zeros(3, 3);
ee = 0.826;
x = true;  % flag, show the direction of axial strain
tolerance = 1e-3;
nstep = 0;
dn = 1e-4;
qLim = 16;
direction = 1;
TestFlag = 1;
LiqFlag = 0;
alphaLiq = zeros(3, 3);
gamma0_max = 0.13;
gamma0 = 0.0;
gamma_bar = zeros(3, 3);
epsv_mono = 0;
gamma_re = [];

while nstep < 50000
    [depsi,direction] = getde(sigma{end}, qLim, dn, direction,TestFlag);
    if nstep == 0
        [depsi,direction] = getde(sigma{end}, qLim, dn / 100, direction,TestFlag);
    end
    G = getG(p(end), ee, G0, pAtmos);
    K = getK(G, nu);
    deE = depsi;
    ds = 2 * G * deE;
    dz = zeros(3, 3);
    dp = 0;
    dAlpha = zeros(3, 3);
    L = 0;

    f = getF(s{end} + ds, alpha, p(end), m);
    r = s{end} / p(end);
    n = getN(r, alpha, m, depsi, LiqFlag);
    cos3Theta = getCos3theta(n);
    g = getg(cos3Theta, c);
    alphaThetaD = getAlphaThetaD(n, p(end), ee, g, M, m, nd, pAtmos, lambdaC, e0, xi, 0, 1);
    alphaThetaB = getAlphaThetaB(n, p(end), ee, g, M, m, nb, pAtmos, lambdaC, e0, xi, 0, 1);
    Ad = getAd(z{end}, n, A0);
    B = getB(cos3Theta, g, c);
    C_val = getC(g, c);
    D = getD(n, Ad, alpha, alphaThetaD);
    RAp = getRAp(B, C_val, n);
    h = getH(ee, p(end), alpha, alphaInit, n, G0, h0, ch, pAtmos);
    Kp = 2 / 3 * p(end) * h * trace((alphaThetaB - alpha) * n);

    if f <= 0
        alphaInit = alpha;
        gamma_re(end+1) = e{end}(1,3);
    elseif p(end) < pMin %&& D > 0
        epsv_mono = epsv_mono + trace(depsi);
        gamma_bar = e{end}+depsi - trace(depsi) / 3 * eye(3);
        gamma_gen = sqrt(2/3 * trace(gamma_bar * gamma_bar));     
        gg = getg(getCos3theta(gamma_bar / sqrt(trace(gamma_bar *gamma_bar))), c);
        dp = 0; 
        if LiqFlag == 0
            %alphaLiq = alpha;            
            LiqFlag = 1;            
            gamma_ex = gamma_gen - gamma0 / gg;
            gamma0 = gamma0 + 0.2 * relu(gamma0_max - gamma0);
            epsv_mono = 0;
        end             
        if trace((depsi - trace(depsi) / 3 * eye(3)) * e{end}) < 0
            gamma_gen = gamma_gen - gamma_ex;
        end
        if gamma_gen * gg + epsv_mono * 1 > gamma0
            dAlpha = alphaThetaD*1.001 - alpha;            
            dp = pMin - p(end);
        end
        ds = dAlpha * pMin;        
        dz = getdz(0, n, z{end}, cz, zmax);
    else
        LiqFlag = 0;
        for i = 1:5
            G = getG(p(end) + dp, ee, G0, pAtmos);
            K = getK(G, nu);
            f = getF(s{end} + ds, alpha + dAlpha, p(end) + dp, m);

            if abs(f) < tolerance
                break;
            end

            % 计算dL
            tmp = -2 * G * RAp + D * K * alpha - 2.0/3 * p(end) * h * (alphaThetaB - alpha);
            s_temp = s{end} + ds - (p(end) + dp) * (alpha + dAlpha);
            denominator = trace(tmp' * s_temp) / sqrt(trace(s_temp' * s_temp)) + sqrt(2/3) * D * K * m;

            dL = -f / denominator;
            L = L + dL;
            
            dp = K * (-relu(L) * D);
            ds = (depsi - relu(L) * RAp) * 2 * G;
            dAlpha = getdAlpha(L, h, alphaThetaB, alpha);            
        end
        dz = getdz(relu(L) * D, n, z{end}, cz, zmax);
    end

    p(end+1) = p(end) + dp;
    s{end+1} = s{end} + ds;
    z{end+1} = z{end} + dz;
    e{end+1} = e{end} + depsi;
    sigma{end+1} = s{end} + p(end) * eye(3);
    alpha = alpha + dAlpha;

    nstep = nstep + 1;

end

e11 = getq(e,TestFlag) * 2;
q = getq(sigma,TestFlag);
p = p';
gamma_re_c = [gamma_re(2)];
for i = 3:length(gamma_re)
    if gamma_re(i) * gamma_re_c(end) < 0
        gamma_re_c(end + 1) = gamma_re(i);
    end
end
gamma_re_c = gamma_re_c';

figure('Position', [200, 200, 1200, 500]); % [left, bottom, width, height]
subplot(1, 2, 1);
plot(p, q, 'r-', 'LineWidth', 1);
xlabel('平均应力 p');
ylabel('偏应力 q');
title('p-q关系图');
grid on;

subplot(1, 2, 2);
plot(e11, q, 'b-', 'LineWidth', 1);
xlabel('轴向应变 e_{11}');
ylabel('偏应力 q');
title('e_{11}-q关系图');
grid on;

function [de,direction] = getde(sigma, qLim, dn, direction, testFlag)

    if testFlag == 1
        q = sigma(1, 1) - sigma(3, 3);
        if q > qLim
            direction = -1;
            dn = dn / 100;
        elseif q < -qLim
            direction = 1;
            dn = dn / 100;
        end
        de = diag([1,-0.5,-0.5]) * dn * direction; 
    elseif testFlag == 2
        q = sigma(1, 3);
        if q > qLim
            direction = -1;
            dn = dn / 100;
        elseif q < -qLim
            direction = 1;
            dn = dn / 100;
        end
        de = [0,0,1;0,0,0;1,0,0] * dn * direction; 
    end
end

function q_val = getq(sigma, testFlag)
    q_val = zeros(length(sigma), 1);
    for i = 1:length(sigma)
        if testFlag == 1
            q_val(i) = sigma{i}(1, 1) - sigma{i}(3, 3);
        elseif testFlag == 2
            q_val(i) = sigma{i}(1, 3);
        end
    end
end

function y = relu(x)
    % ReLU函数
    if x > 0
        y = x;
    else
        y = 0;
    end
end

function G = getG(p, e, G0, pAtmos)
    % 剪切模量
    if p < 0
        fprintf('p: %f\n', p);
        p = 0;
    end
    G = G0 * pAtmos * (2.97 - e)^2 / (1 + e) * sqrt(p / pAtmos);
end

function K = getK(G, nu)
    % 体积模量
    K = 2.0 / 3.0 * (1 + nu) / (1 - 2 * nu) * G;
end

function ec = getEc(pc, lambdaC, e0, xi, pAtmos, eA, A)
    % 临界状态砂土孔隙比
    if pc < 0
        fprintf('pc: %f\n', pc);
        pc = 1e-6;
    end
    ec = e0 - lambdaC * (pc / pAtmos)^xi - eA* (A - 1);    
end

function f = getF(s, alpha, p, m)
    % 是否进入塑性的判别式，小于0为弹性
    temp = s - p * alpha;
    f = sqrt(trace(temp' * temp)) - sqrt(2 / 3.0) * p * m;
end

function n = getN(r, alpha, m, de, LiqFlag)
    % n, tr(n)=0，tr(n^2)=1
    n = (r - alpha) / (sqrt(2 / 3.0) * m);
    if abs(trace(n' * n) - 1) > 1e-5 && LiqFlag == 0
        n_1 = n / sqrt(trace(n' * n));
        if trace(de' * n_1) < 0
            n = -1 * n_1;
        else
            n = n_1;
        end
    elseif abs(trace(n' * n) - 1) > 1e-5 && LiqFlag == 1
        n_1 = de / sqrt(trace(de' * de));
        if trace(de' * n_1) < 0
            n = -1 * n_1;
        else
            n = n_1;
        end
    end
end

function cos3theta = getCos3theta(n)
    % cos3θ
    cos3theta = sqrt(6) * trace(n * n * n);
    if cos3theta > 1.0
        cos3theta = 1.0;
    elseif cos3theta < -1.0
        cos3theta = -1.0;
    end
end

function g = getg(cos3theta, c)
    % g函数
    g = 2 * c / (1 + c - (1 - c) * cos3theta);
end

function alphaThetaD = getAlphaThetaD(n, p, e, g, M, m, nd, pAtmos, lambdaC, e0, xi, eA, A)
    % AlphaThetaD函数
    if p < 0
        fprintf('p: %f\n', p);
        p = 1e-6;
    end
    alphaThetaD = sqrt(2.0 / 3) * (g * M * exp(nd * (e - getEc(p, lambdaC, e0, xi, pAtmos, eA, A))) - m) * n;
end

function B = getB(cos3theta, g, c)
    % B函数
    B = 1 + 3.0 / 2 * (1 - c) / c * g * cos3theta;
end

function C_val = getC(g, c)
    % C函数
    C_val = 3 * sqrt(3.0 / 2) * (1 - c) / c * g;
end

function Ad = getAd(z, n, A0)
    % Ad函数
    Ad = A0 * (1 + relu(trace(z' * n)));
end

function D = getD(n, Ad, alpha, alphaThetaD)
    % D函数
    D = Ad * trace((alphaThetaD - alpha)' * n);
end

function alphaThetaB = getAlphaThetaB(n, p, e, g, M, m, nb, pAtmos, lambdaC, e0, xi, eA, A)
    % AlphaThetaB函数
    if p < 0
        fprintf('p: %f\n', p);
        p = 1e-6;
    end
    alphaThetaB = sqrt(2.0 / 3) * (g * M * exp(-nb * (e - getEc(p, lambdaC, e0, xi, pAtmos, eA, A))) - m) * n;
end

function h = getH(e, p, alpha, alphaInit, n, G0, h0, ch, pAtmos)
    % h函数
    if p < 0
        fprintf('p: %f\n', p);
        p = 1e-6;
    end
    b0 = G0 * h0 * (1 - ch * e) * (p / pAtmos)^(-0.5);
    tmp = trace((alpha - alphaInit)' * n);
    if tmp == 0
        tmp = 1e-8;
    end
    h = b0 / tmp;
end

function RAp = getRAp(B, C_val, n)
    % R'函数
    %RAp = B * n - C_val * (n * n - 1.0 / 3 * eye(3));
    RAp = B * n;
end

function dAlpha = getdAlpha(L, h, alphaThetaB, alpha)
    % dα函数
    if L > 0
        dAlpha = L * 2 / 3 * h * (alphaThetaB - alpha);
    else
        dAlpha = zeros(size(alpha));
    end
end

function dz = getdz(depsvp, n, z, cz, zmax)
    % dz函数
    dz = -cz * relu(-depsvp) * (zmax * n + z);

end

function dF = getdF(depsvp, n, z, cF)
    % dF函数
    dF = cF * depsvp * (n - z);
end

function [ds1, ds2] = getAnswer(s, p, alpha, m)
    % 解方程
    temp = s - p * alpha;
    a = 3 / 2;
    b = 2 * (temp(1,1) - temp(2, 2) / 2 - temp(3, 3) / 2);
    c = temp(1, 1)^2 + temp(2, 2)^2 + temp(3, 3)^2 - 2 / 3 * p^2 * m^2;
    tmp1 = (-b - sqrt(b^2 - 4 * a * c)) / (2 * a);
    tmp2 = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a);
    ds1 = diag([tmp1, -tmp1/2, -tmp1/2]);
    ds2 = diag([tmp2, -tmp2/2, -tmp2/2]);
end

function alpha_val = getAlpha(s, p, flag, m)
    % 直接获得alpha
    a = 3 / 2;
    b = -2 * (s(1,1) - s(2, 2) / 2 - s(3, 3) / 2);
    c = s(1, 1)^2 + s(2, 2)^2 + s(3, 3)^2 - 2 / 3 * p * p * m * m;
    tmp = 4 * p * p * m * m;
    tmp1 = (-b - sqrt(tmp)) / (2 * a * p);
    tmp2 = (-b + sqrt(tmp)) / (2 * a * p);
    if flag == 1
        alpha_val = diag([tmp1, -tmp1/2, -tmp1/2]);
    else
        alpha_val = diag([tmp2, -tmp2/2, -tmp2/2]);
    end
end

function [ds, L] = getDrainedds(Kp, n, r, G, RAp, de)
    % 直接求解ds
    dd = de(1, 1);
    if Kp == 0
        Kp = 1e-7;
    end
    dp = dd / (1 / G + RAp(1,1) / Kp * (n(1, 1) * 2 - n(2, 2) - n(3, 3) - trace(n' * r)));
    ds = diag([2 * dp, -dp, -dp]);
    L = (trace(n' * ds) - trace(n' * r) * dp) / Kp;
    if L < 0
        ds = de * 2 * G;
        dp = ds(1, 1) / 2;
        L = (trace(n' * ds) - trace(n' * r) * dp) / Kp;
    end
end