function [rho, phi, segments, value] = findSymmetry(inputImage, parameters)

%%
if nargin < 2
    parameters = get_default_parameters()
end

img = double(rgb2gray(inputImage));

verbose = parameters.verbose;
numberOfLines = parameters.numberOfLines;
symmetryAngles = parameters.symmetryAngles;

%% Scan angles
lineSet = cell(numel(symmetryAngles),1);

if verbose > 0; disp('Start symmetry scanning'); end;
if verbose > 0; tic; end;

% If the Parallel computing toolbox is implemented, the parfor will run
% using a pool of worker threads. It should run like an ordinary for-loop
% otherwize.
parfor phiIndx = 1:1:numel(symmetryAngles)
    phi = symmetryAngles(phiIndx);
    if verbose > 0; fprintf('Phi: %.2fpi\n',phi/pi); end;
    
    [SYM, houghRhoOffset, lengthOffset] = computeSymmetryMatrix(img, phi, parameters);
        
    [rhos, values, lowerBounds, upperBounds] = selectCandidate(SYM, parameters);

    rhos = rhos + houghRhoOffset;
    lowerBounds = lowerBounds + lengthOffset;
    upperBounds = upperBounds + lengthOffset;

    lines = [rhos; values; lowerBounds; upperBounds];                            
    lineSet{phiIndx} = [phi * ones(1,size(lines,2)); lines];
        
end
if verbose > 0; toc; end;

allLines = cat(2,lineSet{:});

[~, i] = sort(allLines(3,:));
i = i(end:-1:1);
allLinesSorted = allLines(:,i);

numberOfLines = min(size(allLinesSorted,2), numberOfLines);

phi = allLinesSorted(1,1:numberOfLines);
rho = allLinesSorted(2,1:numberOfLines);
lo = allLinesSorted(4,1:numberOfLines);
hi  = allLinesSorted(5,1:numberOfLines);
value  = allLinesSorted(3,1:numberOfLines);


%% Compute segments
segments = cell(numberOfLines,1);
for i = 1:numberOfLines
    theta = phi(i) - pi/2;
    [cx, cy] = pol2cart(theta,rho(i));
    [lx, ly] = pol2cart(theta+pi/2,lo(i));
    [hx, hy] = pol2cart(theta+pi/2,hi(i));

    segments{i} = [[cx+lx;cx+hx], [cy+ly;cy+hy]];
end

end






