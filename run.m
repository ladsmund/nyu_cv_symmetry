clearvars
close all
%%
addpath('symmetryDectection')

root = './SymmetryDBppGT/';
inputFolders = {'S/','M/'};
outputRoot = './SymmetryDBppGT/symmetry/';
% [~,~] = mkdir(outputFolder);

%%
% symmetryAngles = ((0:16)/16)*pi;
symmetryAngles = ((0:8)/8)*pi;
searchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];
% searchAngles = [[-pi/4;pi/4], [0;0], [pi/4;-pi/4]];

% ranges = [2:2:100];
ranges = [2:4:50];
numberOfLines = 10;
sigmas = [2];
% sigmas = [2 4 10];

%%
for sigma = sigmas
    fprintf('Sigma = %.2f\n',sigma);
    tic
    for inputFolder = inputFolders
        folder = [root inputFolder{1}];
        outputFolder = [outputRoot num2str(sigma) inputFolder{1}];
        [~,~] = mkdir(outputFolder);

        imgType = '*.png';
        images = dir([folder imgType]);
        
        for idx = 1:length(images)
            
            path = [folder images(idx).name];
            [~, name, ext] = fileparts(path);        
            disp(path)

            img = imread(path);

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
                
                %% Save data file with symmetry information
                save(sprintf('%s.mat',[outputFolder name]), ...
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
                     'sigmas');

                    %% Create an image file with the symmetry lines
                    fig = figure(1);
                    fig.Visible = 'off';
                    imshow(img)
                    colors = hot(3);
                    hold on
                    colorIndx = 1;
                    for i = 1:min(numberOfLines,numel(phi))
                        theta = phi(i) - pi/2;

                        [cx, cy] = pol2cart(theta,rho(i));

                        [lx, ly] = pol2cart(theta+pi/2,lo(i));
                        [hx, hy] = pol2cart(theta+pi/2,hi(i));

                        [dx, dy] = pol2cart(theta+pi/2,2000);
                        line([cx+lx,cx+hx], [cy+ly,cy+hy],'Color',colors(colorIndx,:), 'LineWidth', 4)
                        colorIndx = colorIndx + 1;
                    end
                    hold off

                    frame = getframe(fig);
                    imwrite(frame.cdata,[outputFolder, name, ext])

            catch ME            
                disp('Error');  
                disp(ME.message);
                toc
                continue;
            end
            toc    
        end        
    end
end