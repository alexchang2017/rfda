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
sprintf('%.7f, ', interp1(x, [y; y]', xi))
sprintf('%.7f, ', interp1(x, [y; y]', xi, 'spline'))

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

save('testMat/interp2_res.mat', 'V_1_l', 'V_1_s', 'V_2_l', 'V_2_s', 'V_3_l', 'V_3_s')

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

% locLinearRotate2d
x = 0:0.1:1;
[x1, x2] = meshgrid(x, x);
x = [x1(:), x2(:)];
idx = x(:,1) > x(:,2);
x = x(x(:,1) ~= x(:,2), :);
y = x .^ 2 .* 2 + 3 .* x + 3;
y = y(:,1) .* (y(:,2) .* 0.5 - 1);
outMat = [0:0.2:1; 0:0.2:1];
cnt2 = [ones(1, 55), ones(1, 55) * 2];
w = ones(1,length(y));
w2 = [ones(1, 55), ones(1, 55) * 2];
[~, est] = rotate_mlwls([1, 1], 'gauss', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([2, 2], 'gauss', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'gausvar', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([2, 2], 'gausvar', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'epan', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([2, 2], 'epan', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'quar', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([2, 2], 'quar', x', y, w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'gauss', x', y ./ cnt2', w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'gauss', x', y, w2, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'epan', x', y ./ cnt2', w, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([1, 1], 'epan', x', y, w2, outMat, 1);
sprintf('%.6f, ', est)
[~, est] = rotate_mlwls([0.075, 0.075], 'epan', x(idx, :)', y(idx), w(idx), outMat, 1);
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
% bwCandChooser2
h0 = getMinb(t, out1, regular) * 0.2;
q = (r/(4*h0))^(1/9);
bwc = sort(q.^(0:9).*h0);
sprintf('%.6f, ', bwc) 

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
% bwCandChooser2
h0 = getMinb(t, out1, regular) * 0.2;
q = (r/(4*h0))^(1/9);
bwc = sort(q.^(0:9).*h0);
sprintf('%.6f, ', bwc) 

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
% bwCandChooser2
h0 = getMinb(t, out1, regular) * 0.2;
q = (r/(4*h0))^(1/9);
bwc = sort(q.^(0:9).*h0);
sprintf('%.6f, ', bwc) 

%% test results for rawCov_res (not used)
for example_case = 0:2
  switch example_case
    case 0
      load('exData/sparseExData.mat')
      exampleName = 'sparse.mat';
    case 1
      load('exData/irregularExData.mat')
      exampleName = 'irregular.mat';
    case 2
      load('exData/regularExData.mat')
      exampleName = 'regular.mat';
  end

  tt = cell2mat(t);
  yy = cell2mat(y);
  w = ones(1, length(yy));
  bopt = gcv_lwls(yy, tt, 'gauss', 1, 1, 0, example_case, 'on');
  bopt = sqrt(minb(tt,2)*adjustBW1('gauss',bopt, 1, 0, example_case,'on'));
  [~, muDense] = lwls(bopt, 'gauss', 1, 1, 0, tt, yy', w, sort(unique(tt)), 0);
  rawCov1 = getRawCov(y, t, sort(unique(tt)), muDense, example_case, 0);
  rawCov2 = getRawCov(y, t, sort(unique(tt)), muDense, example_case, 1);
  % save(['testMat/rawCov_res_', exampleName], 'rawCov1', 'rawCov2') 
end

%% test for getRawCov and locLinear2d / locLinearRotate2d
ngrid = 30;
for example_case = 0:2
  switch example_case
  case 0
    regular = 0;
    y = {[0.1,0.2,0.3], [0.4,0.5,0.6], [0.7,0.8], [0.9,1]};
    t = {[1,2,4], [3,5,7], [6,10], [8,9]};
  case 1
    regular = 1;
    y = {[0.1,0.2,0.3,0.4], [0.5,0.6,0.7], [0.8,0.9,1], [1.1,1.2,1.3]};
    t = {[1,2,3,4], [1,2,4], [2,3,4], [1,3,4]};
  case 2
    regular = 2;
    y = {[0.1,0.2,0.3], [0.4,0.5,0.6], [0.7,0.8,0.9], [0.9,1,1]};
    t = {[1,2,3], [1,2,3], [1,2,3], [1,2,3]};
  end
  
  tt = cell2mat(t);
  yy = cell2mat(y);
  out1 = sort(unique(tt));
  bopt = gcv_lwls(yy, tt, 'gauss', 1, 1, 0, regular, 'on');
  [~, mu] = lwls(adjustBW1('gauss',bopt, 1, 0, 0,'on'), 'gauss', 1, 1, 0, tt, yy', ones(1,length(tt)), out1, 0);
  rcov = getRawCov(y, t, out1, mu, regular, 0);
  gcv_bw = gcv_mullwlsn(t, 30, regular, 1, 'gauss', rcov, 'on');

  out21 = linspace(min(tt),max(tt),ngrid);
  tneq = find(rcov.tpairn(1,:) ~= rcov.tpairn(2,:));
  win = ones(1,length(tneq));
  if example_case ~= 1
    [~, xcov_10]= mullwlsk([1, 1], 'gauss', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21);
    [~, xcov_15]= mullwlsk([1.5, 1.5], 'gauss', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21);
    [~, xcov_20]= mullwlsk([2, 2], 'gauss', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21);
  else
    [~, xcov_10]= mullwlsk([1, 1], 'gauss', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21, rcov.count(tneq));
    [~, xcov_15]= mullwlsk([1.5, 1.5], 'gauss', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21, rcov.count(tneq));
    [~, xcov_20]= mullwlsk([2, 2], 'gauss', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21, rcov.count(tneq));
  end
  
  eval(['rcov_case', int2str(regular), '=[rcov.tpairn; rcov.cxxn]'';'])
  eval(['xcov_10_case', int2str(regular), '=xcov_10;'])
  eval(['xcov_15_case', int2str(regular), '=xcov_15;'])
  eval(['xcov_20_case', int2str(regular), '=xcov_20;'])
  eval(['gcv_bw_case', int2str(regular), '=gcv_bw;'])
  
  if example_case == 2
    [~, xcov_30_gaussvar] = mullwlsk([3, 3], 'gausvar', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21);
    [~, xcov_30_epan] = mullwlsk([3, 3], 'epan', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21);
    [~, xcov_30_quar] = mullwlsk([3, 3], 'quar', rcov.tpairn(:,tneq), rcov.cyy(tneq)', win, out21, out21);
  end
end
save('testMat/covRes.mat', 'rcov_case0', 'rcov_case1', 'rcov_case2', 'xcov_10_case0', ...
  'xcov_10_case1', 'xcov_10_case2', 'xcov_15_case0', 'xcov_15_case1', 'xcov_15_case2', ...
  'xcov_20_case0', 'xcov_20_case1', 'xcov_20_case2', 'xcov_30_gaussvar', 'xcov_30_epan', ...
  'xcov_30_quar', 'gcv_bw_case0', 'gcv_bw_case1', 'gcv_bw_case2');

