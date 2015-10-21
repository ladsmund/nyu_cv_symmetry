clearvars
close all
%%
addpath('symmetryDectection')
%% Read image files

% image = 'I120'; % Glass
image = 'I174'; % Helicopter
% image = 'I172'; % Plane difficult
% image = 'I175'; % Helicopter

folder = '~/workspace/SymmetryDBpp/S/';
image_path = sprintf('%s%s.png',folder,image);
img = imread(image_path);


%%
parameters = get_default_parameters()

% symmetryAngles = ((1:16)/16)*pi;
parameters.symmetryAngles = ((1:4)/4)*pi;
% symmetryAngles = pi/2;
% searchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];
% searchAngles = [[-pi/4;pi/4], [0;0], [pi/4;-pi/4]];
% searchAngles = [-pi/4,0, pi/4];
searchAngles = [-pi/4, 0,pi/4];
parameters.searchAngles = [searchAngles - pi/2; pi/2- searchAngles];


%
% searchAngles = [[pi/2;pi/2]];
% searchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];

parameters.numberOfLines = 10;

parameters.sigmas = 2;

% ranges = [2:2:50];
parameters.ranges = [2:2:50];


%%
[rho, phi, segments] = findSymmetry(img, parameters);

%

data = struct();
data.segments = segments;
data.image_path = image_path;

figure(1);
show_result(data)

