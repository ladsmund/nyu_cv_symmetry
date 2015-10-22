function [rho, phi, segments, value] = findSymmetry(inputImage, parameters)

%%

if nargin < 2
    parameters = get_default_parameters()
end


img = double(rgb2gray(inputImage));
imgHeight= size(img,1);
imgWidth = size(img,2);

verbose = parameters.verbose;
visualize = parameters.visualize;
searchRange = parameters.searchRange;
numberOfLines = parameters.numberOfLines;
symmetryAngles = parameters.symmetryAngles;
filterCombinator = parameters.filterCombinator;
symmetryMetric = parameters.symmetryMetric;
sigmas = parameters.sigmas;

searchAngles = parameters.searchAngles;
searchAngles = [searchAngles - pi/2; pi/2- searchAngles];

%%

biasFactor = parameters.distanceBiasAlpha * searchRange.^-2;
biasFactorMatrix = repmat(permute(biasFactor,[3,1,2]),imgHeight,imgWidth,1);

%%

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
    
    %% Compute Symmetry matrix
    % The similarity matrix is a 3D matrix where the coordinates are the with,
    % height and search distance.
    SIM = zeros(imgHeight,imgWidth,length(searchRange));
    for i = 1:size(searchAngles,2)
        theta = searchAngles(:,i);    
        theta1 = mod(phi+theta(1),2*pi);
        theta2 = mod(phi+theta(2),2*pi);
        
        J1 = morletFilter(img,theta1, sigmas);
        J2 = morletFilter(img,theta2, sigmas);
        
        for rIndx = 1:numel(searchRange)
            r = searchRange(rIndx);
            [dx, dy] = pol2cart(phi-pi/2,r);
            d = ([dx, dy]);
            J1t = imtranslate(J1,d);
            J2t = imtranslate(J2,-d);
            
            % Combining the two filtered images
            SIM(:,:,rIndx) = SIM(:,:,rIndx) + filterCombinator(J1t, J2t);
        end                
    end

    % Compute real values symmetry metric for voting    
    SIM = symmetryMetric(SIM);
    SIM = SIM + biasFactorMatrix;

    %% Transform into rho/line/distance space
    % By rotating the similarity matrix with respect to the houghAngle.
    % The function imrotate used degree and rotates in the opposite direction.
    houghAngle = phi - pi/2;
    rotAngle= 180*(houghAngle)/pi;
    imgRot = imrotate(SIM,rotAngle,'bilinear');

    houghRhoOffset = min(0,round(imgHeight * sin(houghAngle)));

    lenOffset = min(0,-round(imgWidth * sin(houghAngle)));
    lengthOffset = 0;

    [rhos, values, lowerBounds, upperBounds] = selectCandidate(imgRot, parameters);

    rhos = rhos + houghRhoOffset;
    lowerBounds = lowerBounds + lengthOffset + lenOffset;
    upperBounds = upperBounds + lengthOffset + lenOffset;

    lines = [rhos; values; lowerBounds; upperBounds];                            
    lineSet{phiIndx} = [phi * ones(1,size(lines,2)); lines];
        
end
if verbose > 0; toc; end;
%% Generate matrix for hough column

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






