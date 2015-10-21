function output_files = batch_find_symmetry(input_path, output_path, parameters)
%%
addpath('symmetryDectection')
if nargin < 3
    parameters = get_default_parameters();
end

%%
fprintf('Starting batch find symmetry\n')
fprintf('input_path: %s\n', input_path)
fprintf('output_path: %s\n', output_path)
disp(parameters)

%%

images = dir(input_path);
folder = fileparts(input_path);
image_number = 0;
sigma = parameters.sigmas;
image_number = image_number + 1;

output_files = cell(numel(images), 1);
tic

% outputFolder = [output_path sprintf('sigma_%s/',num2str(sigma))];
outputFolder = [output_path];
[~,~] = mkdir(outputFolder);

for idx = 1:length(images)
    image_path = [folder '/' images(idx).name];
    [~, name, ext] = fileparts(image_path); 

    img = imread(image_path);

    try
        [rho, phi, segments] = findSymmetry(img, parameters);            

        %% Save data file with symmetry information
        image_output_path = sprintf('%s.mat',[outputFolder name]);
        output_files{image_number} = image_output_path;

        save(image_output_path, ...
             'parameters',...
             'image_path',...
             'folder', ...
             'name', ...
             'ext', ...
             'rho', ...
             'phi', ...
             'segments');

    catch ME            
        disp('Error');  
        disp(ME.message);
        toc
        continue;
    end    
    toc
end

%%
end