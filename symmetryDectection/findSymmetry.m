function [rho, phi, segments] = findSymmetry(inputImage, varargin)
%% Constants
FLOAT_EQUALITY_PRECITION = 8;

if verLessThan('matlab', '8.0.0')
    % I don't know when matlab introduced precition for round
    round_angle = @(x) round(x * 10^FLOAT_EQUALITY_PRECITION) / 10^FLOAT_EQUALITY_PRECITION;
else
    round_angle = @(x) round(x,FLOAT_EQUALITY_PRECITION);    
end

%% Parse input
p = inputParser;

defaultSigmas = 4;
defaultNumberOfSymmetryAngles = 32;
defaultSearchRange = 10:5:50;
defaultSearchAngles = [[-pi/4;pi/4], [pi/4;-pi/4]];
% defaultSearchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];
defaultNumberOrLines = 5;
defaultVerbose = 1;
defaultVisualize = 0;

addRequired(p,'image',@isnumeric);
addOptional(p,'sigmas',defaultSigmas,@isnumeric)
addOptional(p,'numberOfSymmetryAngles',defaultNumberOfSymmetryAngles,@isnumeric)
addOptional(p,'searchRange',defaultSearchRange,@isnumeric)
addOptional(p,'searchAngles',defaultSearchAngles,@isnumeric)
addOptional(p,'numberOfLines',defaultNumberOrLines,@isnumeric)
addOptional(p,'verbose',defaultVerbose,@isnumeric)
addOptional(p,'visualize',defaultVisualize,@isnumeric)
addOptional(p,'symmetryAngles',null(0), @isnumeric)

parse(p,inputImage,varargin{:});
args = p.Results;

img = double(rgb2gray(args.image));
imgHeight= size(img,1);
imgWidth = size(img,2);

if ~numel(args.symmetryAngles)
    symmetryAngles = (0:args.numberOfSymmetryAngles-1) * (2*pi / args.numberOfSymmetryAngles);
else
    symmetryAngles = args.symmetryAngles;
end

verbose = args.verbose;
visualize = args.visualize;
searchRange = args.searchRange;
numberOfLines = args.numberOfLines;
searchAngles = args.searchAngles;
filterAngles = [];
%%
for sa = searchAngles(:)'
    filterAngles = [filterAngles mod(symmetryAngles+sa,2*pi)];
    filterAngles = unique(round_angle(filterAngles));
end
filterAngles = sort(squeeze(filterAngles));

%% Pregenerate a set of filtered images
if verbose > 0; disp('Generate set of filtered images');end;
if verbose > 0; tic; end;
filterSet = generateFilterSet(img, filterAngles, args.sigmas,visualize);
if verbose > 0; toc; end;

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
      
    imageIndexes = zeros(2,size(searchAngles,2));
    for i = 1:size(searchAngles,2)
        theta = searchAngles(:,i);    
        theta1 = round_angle(mod(phi+theta(1),2*pi));
        theta2 = round_angle(mod(phi+theta(2),2*pi));        
        imageIndexes(:,i) = [find(filterAngles == theta1); ...
                             find(filterAngles == theta2)];
    end
    
    %% Compute similarity matrix
    % The similarity matrix is a 3D matrix where the coordinates are the with,
    % height and search distance.
    SIM = zeros(imgHeight,imgWidth,length(searchRange));
    for i = imageIndexes
        tIndx1 = i(1);
        tIndx2 = i(2);

        J1 = filterSet.filtered{tIndx1};
        J2 = filterSet.filtered{tIndx2};

        for rIndx = 1:numel(searchRange)
            r = searchRange(rIndx);
            [dx, dy] = pol2cart(phi-pi/2,r);
            d = ([dx, dy]);
            J1t = imtranslate(J1,d);
            J2t = imtranslate(J2,-d);
            
            % Combining the two filtered images
            SIM(:,:,rIndx) = SIM(:,:,rIndx) + (J1t .* conj(J2t));
        end                
    end

    % Compute real values symmetry metric for voting
    SIM = real(sqrt(SIM .* conj(SIM)));

    %% Transform into rho/line/distance space
    % By rotating the similarity matrix with respect to the houghAngle.
    % The function imrotate used degree and rotates in the opposite direction.
    houghAngle = phi - pi/2;
    rotAngle= 180*(houghAngle)/pi;
    imgRot = imrotate(SIM,rotAngle,'bilinear');

    houghRhoOffset = min(0,round_angle(imgHeight * sin(houghAngle)));

    lenOffset = min(0,-round_angle(imgWidth * sin(houghAngle)));
    lengthOffset = 0;

    [rhos, values, lowerBounds, upperBounds] = selectCandidate(imgRot,numberOfLines,args.sigmas,0);

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






