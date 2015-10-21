clearvars
close all
%%

% image = 'I120'; % Glass
% image = 'I174'; % Helicopter
image = 'I172'; % Plane

input_path = '~/workspace/SymmetryDBpp/S/';
output_path = '~/workspace/SymmetryDBpp/output/S/';
gt_path = '~/workspace/SymmetryDBpp/S/';

%%
parameters = get_default_parameters();
parameters.symmetryAngles = ((0:16)/16)*pi;
% parameters.symmetryAngles = ((0:8)/8)*pi;
% parameters.symmetryAngles = pi/2;
parameters.searchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];
% parameters.searchAngles = [[-pi/4;pi/4], [0;0], [pi/4;-pi/4]];

% parameters.ranges = [2:2:100];
parameters.ranges = [2:4:80];
% parameters.ranges = [2:4:50];
parameters.numberOfLines = 20;
parameters.sigmas = [2];
% parameters.sigmas = [2 4 10];

%%
output_files = batch_find_symmetry([input_path image '.png'], output_path, parameters);

%%

figure(1)

% image = 'I171';
output_path = '~/workspace/SymmetryDBpp/output/S/sigma_2/';

subplot(1,2,1)
show_result(load([gt_path image '.mat']), [gt_path image '.png']);

subplot(1,2,2)
% show_result(output_files{1});
show_result(load([output_path image '.mat']));
