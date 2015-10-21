clearvars
close all
%%
addpath('symmetryDectection')
%%
image = 'I120'; % Glass
% image = 'I174'; % Helicopter
% image = 'I172'; % Plane

input_path = '~/workspace/SymmetryDBpp/S/';
output_path = '~/workspace/SymmetryDBpp/output/S/sigma2/';
gt_path = '~/workspace/SymmetryDBpp/S/';

%%
parameters = get_default_parameters();
parameters.symmetryAngles = ((1:4)/4)*pi;
% parameters.symmetryAngles = pi/2;
parameters.searchAngles = [-pi/4, 0, pi/4];

parameters.ranges = [2:4:80];

parameters.numberOfLines = 10;
parameters.sigmas = [2];
parameters.verbose = 1;

%%
output_files = batch_find_symmetry([input_path image '.png'], output_path, parameters);

%%

figure(1)

% image = 'I171';
% output_path = '~/workspace/SymmetryDBpp/output/S/sigma2/';

subplot(1,2,1)
show_result(load([gt_path image '.mat']), [gt_path image '.png']);

subplot(1,2,2)
show_result(load(output_files{1}));
% show_result(load([output_path image '.mat']));
