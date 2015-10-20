function output_files = batch_find_symmetry(input_path, output_path)
%%

% input_path = '~/workspace/SymmetryDBpp/S/I096.png';
% output_path = '~/workspace/SymmetryDBpp/output/S/';

%%
addpath('symmetryDectection')

% symmetryAngles = ((0:16)/16)*pi;
symmetryAngles = ((0:8)/8)*pi;
searchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];
% searchAngles = [[-pi/4;pi/4], [0;0], [pi/4;-pi/4]];

% ranges = [2:2:100];
ranges = [2:2:80];
% ranges = [2:4:50];
numberOfLines = 20;
sigmas = [2];
% sigmas = [2 4 10];

%%

output_files = cell(numel(sigmas), 1);

images = dir(input_path);
folder = fileparts(input_path);
image_number = 0;
for sigma = sigmas
    image_number = image_number + 1;
    fprintf('Sigma = %.2f\n',sigma);
    tic
    
    outputFolder = [output_path sprintf('sigma_%s/',num2str(sigma))];
    [~,~] = mkdir(outputFolder);

    for idx = 1:length(images)
        image_path = [folder '/' images(idx).name];
        [~, name, ext] = fileparts(image_path); 

        img = imread(image_path);

        try
            [rho, phi, lo, hi] = findSymmetry(...
                img...
                ,'visualize',0 ...
                ,'verbose',0 ...
                ,'searchRange',ranges ...
                ,'sigmas',sigma ...
                ,'numberOfLines',numberOfLines ...
                ,'searchAngles',searchAngles ...
                ,'symmetryAngles',symmetryAngles ...
                );

            %% Compute segments
            number_of_segments = min(numberOfLines,numel(phi));
            segments = cell(number_of_segments,1);
            for i = 1:number_of_segments
                theta = phi(i) - pi/2;
                [cx, cy] = pol2cart(theta,rho(i));
                [lx, ly] = pol2cart(theta+pi/2,lo(i));
                [hx, hy] = pol2cart(theta+pi/2,hi(i));
                
                segments{i} = [[cx+lx;cx+hx], [cy+ly;cy+hy]];
                
            end
            
            
            %% Save data file with symmetry information
            image_output_path = sprintf('%s.mat',[outputFolder name]);
            output_files{image_number} = image_output_path;
            
            save(image_output_path, ...
                 'image_path',...
                 'folder', ...
                 'name', ...
                 'ext', ...
                 'rho', ...
                 'phi', ...
                 'lo', ...
                 'hi', ...
                 'symmetryAngles', ...
                 'searchAngles', ...
                 'ranges', ...
                 'numberOfLines', ...
                 'segments',...
                 'sigma');

        catch ME            
            disp('Error');  
            disp(ME.message);
            toc
            continue;
        end    
        toc
    end

end

%%
end