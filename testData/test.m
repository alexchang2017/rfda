addpath('../../PACE2.17/PACE')

%% other
% trapz
xi = 1:5;
x = xi .^ 2;
x2 = [xi .^2; xi .^3];
sprintf('%.4f', trapz(x))
sprintf('%.4f', trapz(xi, x))
sprintf('%.4f, ', trapz(x2'))
sprintf('%.4f, ', trapz(xi, x2'))

% interp1
x = [0.8, 0.3, 0.1, 0.6];
y = x .^ 2;
xi = linspace(0.1, 0.8, 5);
sprintf('%.7f, ', interp1(x, y, xi))
sprintf('%.7f, ', interp1(x, y, xi, 'spline'))
sprintf('%.7f, ', interp1([0.1, 0.9], [0.3, 0.8], [0.1, 0.5, 0.9], 'spline'))

x = [0.8, 0.3, 0.1, 0.6, 0.9, 0.5, 0.2, 0.0, 0.7, 1.0, 0.4];
y = x .^ 2;
xi = linspace(0, 1, 21);
sprintf('%.7f, ', interp1(x, y, xi))
sprintf('%.7f, ', interp1(x, y, xi, 'spline'))

x = 0:(pi/4):(2*pi);
y = sin(x);
xi = 0:(pi/16):(2*pi);
sprintf('%.7f, ', interp1(x, y, xi))

% interp2
A = [13,-1,12;5,4,3;1,6,2];
x = [0,1,4];  y = [10,11,12];
xi = linspace(min (x), max (x), 17);
yi = linspace(min (y), max (y), 26);

V_1_l = interp2(x,y,A,xi,yi', 'linear');
V_1_s = interp2(x,y,A,xi,yi', 'spline');

[X,Y] = meshgrid(-3:3);
V = peaks(X,Y);
[Xq,Yq] = meshgrid(-3:0.25:3);
V_2_l = interp2(X,Y,V,Xq,Yq, 'linear');
V_2_s = interp2(X,Y,V,Xq,Yq, 'spline');

x = [0.1, 0.2, 0.8];
y = [0.1, 0.3, 0.6, 0.8];
v = [0.2,0.3,0.4,0.5;0.5,0.4,0.4,0.4;0.2, 0.4, 0.5, 0.4]';
xi = linspace(0.1, 0.8, 8);
yi = linspace(0.1, 0.8, 8);
V_3_l = interp2(x,y,v,xi,yi', 'linear');
V_3_s = interp2(x,y,v,xi,yi', 'spline');

save('testResult/interp2_res.mat', 'V_1_l', 'V_1_s', 'V_2_l', 'V_2_s', 'V_3_l', 'V_3_s')

% locPoly1d
% - same weight
x = 0:0.1:1;
y = x .^ 2 .* 2 + 3 .* x;
w = ones(1,length(x));
[~, est] = lwls(0.25, 'epan', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 
[~, est] = lwls(0.25, 'quar', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 
[~, est] = lwls(0.25, 'gauss', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 
[~, est] = lwls(0.25, 'gausvar', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 
% - different weight
w = [ones(1,length(x)-2), 0, 0.5];
[~, est] = lwls(0.25, 'epan', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est)
[~, est] = lwls(0.25, 'quar', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 
[~, est] = lwls(0.25, 'gauss', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 
[~, est] = lwls(0.25, 'gausvar', 1, 1, 0, x, y', w, x);
sprintf('%.6f, ', est) 

%% regular
load('exData\regularExData.mat')
regular = 2;
tt = cat(2, t{:});
yy = cell2mat(y);
out1 = sort(unique(tt));
r = range(tt);
h0 = minb(tt, 2)*1.5;
q = (r/(4*h0))^(1/9);
bwc = sort(q.^(0:9).*h0);
% bwCandChooser
sprintf('%.6f, ', bwc)
% gcv_locPoly1d
bopt = gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'gausvar', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'epan', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'quar', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt)
% adjGcvBw1d
adj_bopt = adjustBW1('gauss', gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'off'), 1, 0, regular, 'off');
sprintf('%.6f', adj_bopt)
adj_bopt = adjustBW1('epan', gcv_lwls(yy, tt, 'epan', 1, 1, 0, regular, 'off'), 1, 0, regular, 'off');
sprintf('%.6f', adj_bopt)

%% irregular
load('exData\irregularExData.mat')
regular = 1;
tt = cat(2, t{:});
yy = cell2mat(y);
out1 = sort(unique(tt));
r = range(tt);
h0 = minb(tt, 2)*2;
q = (r/(4*h0))^(1/9);
bwc = sort(q.^(0:9).*h0);
% bwCandChooser
sprintf('%.6f, ', bwc)
% gcv_locPoly1d
bopt = gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'gausvar', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'epan', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'quar', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
% adjGcvBw1d
adj_bopt = adjustBW1('gauss', gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'off'), 1, 0, regular, 'off');
sprintf('%.6f', adj_bopt)
adj_bopt = adjustBW1('epan', gcv_lwls(yy, tt, 'epan', 1, 1, 0, regular, 'off'), 1, 0, regular, 'off');
sprintf('%.6f', adj_bopt)

%% sparse
load('exData\sparseExData.mat')
regular = 0;
tt = cat(2, t{:});
yy = cell2mat(y);
out1 = sort(unique(tt));
r = range(tt);
h0 = minb(tt, 3)*2.5;
q = (r/(4*h0))^(1/9);
bwc = sort(q.^(0:9).*h0);
% bwCandChooser
sprintf('%.6f, ', bwc) 
% gcv_locPoly1d
bopt = gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'gausvar', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'epan', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
bopt = gcv_lwls(yy, tt, 'quar', 1, 1, 0, regular, 'off');
sprintf('%.6f', bopt) 
% adjGcvBw1d
adj_bopt = adjustBW1('gauss', gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'off'), 1, 0, regular, 'off');
sprintf('%.6f', adj_bopt)
adj_bopt = adjustBW1('epan', gcv_lwls(yy, tt, 'epan', 1, 1, 0, regular, 'off'), 1, 0, regular, 'off');
sprintf('%.6f', adj_bopt)
