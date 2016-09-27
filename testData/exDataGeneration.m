function [] = exDataGeneration()
exDataGeneration_f(0);
exDataGeneration_f(1);
exDataGeneration_f(2);
end

function [] = exDataGeneration_f(example_case)
% input: example_case:
%   0 for the sparse case
%   1 for the regular cases
%   2 for the regular with missing cases

% random seed
s = RandStream.create('twister', 'seed', 123);
RandStream.setGlobalStream(s);

ncohort = 200; % 200 subjects in the simulated data
lint = 10;     % maximum time
y = cell(1, ncohort);
t = cell(1, ncohort);

switch example_case
  case 0
    % at most 8 repeated measurements in the simulated data for the sparse case
    mtp = 8;
    example_case_name = 'sparse';
  case {1, 2}
    % 20 repeated measurements in the simulated data for
    % regular and regular with missing cases
    mtp = 20;
    if example_case == 1
      example_case_name = 'regular';
    else
      example_case_name = 'irregular';
    end
end

% Case iii) regular data with missing values (regular = 1)
if example_case == 2
  % generate the number of missing data for each subject
  ni_missing = zeros(1,ncohort);
  for i = 1:ncohort
    ni_missing(i) = poissrnd(randsample(floor(mtp/2), 1), 1, 1);
    if ni_missing(i) >= mtp
      ni_missing(i) = mtp-1;
    end
  end
end

mu_true = @(t) t+sin(t);
phi1_true = @(t, lint) -sqrt(2/lint) * cos(2*pi*t/lint);
phi2_true = @(t, lint) sqrt(2/lint) * sin(2*pi*t/lint);
xi = zeros(ncohort, 2);
for i = 1:ncohort
  ntp=ceil(mtp*rand(1));
  switch example_case
    case 0
      %   Case i) Sparse and irregular case (regular = 0)
      t{i} = lint* rand(1,ntp);
    case 1
      %   Case ii) complete balance case (regular = 2)
      t{i} = linspace(0, lint, mtp);
    case 2
      %   Case iii) regular data with missing values (regular = 1)
      t{i} = linspace(0,lint,mtp);
      if ni_missing(i) > 0
        id_missing = randsample(mtp, ni_missing(i));
        t{i}(id_missing) = [];
      end
  end

  % generate 2 Principal components, 1st PC score: N(0,9), 2nd PC score: N(0,4)
  xi(i, :) = [3*randn(1), 2*randn(1)];
  % generate the repeated measurements with measurement errors
  y{i} = mu_true(t{i}) + xi(i,:) * [phi1_true(t{i}, lint); phi2_true(t{i}, lint)] + ...
    randn(1, length(t{i}));
  % true mean function: t+sin(t)
  % 1st eigenfunction: -sqrt(0.2)*cos(2*t*pi/10)
  % 2nd eigenfuncton: sqrt0.2)*sin(2*t*pi/10)
  % measurement error is distributed as N(0,1)
end
save(['exData/', example_case_name, 'ExData.mat'], 'y', 't')
end
