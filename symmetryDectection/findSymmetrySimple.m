function [rho, phi, segments] = findSymmetry(inputImage, varargin)
%% Constants
FLOAT_EQUALITY_PRECITION = 8;

%% Parse input
p = inputParser;

defaultSigmas = 4;
defaultNumberOfSymmetryAngles = 32;
defaultSearchRange = 10:5:50;
defaultSearchAngles = [[-pi/4;pi/4], [pi/2;pi/2], [pi/4;-pi/4]];
defaultNumberOrLines = 5;
defaultVerbose = 1;
defaultVisualize = 0;

addRequired(p,'image',@isnumeric);
addParameter(p,'sigmas',defaultSigmas,@isnumeric)
addParameter(p,'numberOfSymmetryAngles',defaultNumberOfSymmetryAngles,@isnumeric)
addParameter(p,'searchRange',defaultSearchRange,@isnumeric)
addParameter(p,'searchAngles',defaultSearchAngles,@isnumeric)
addParameter(p,'numberOfLines',defaultNumberOrLines,@isnumeric)
addParameter(p,'verbose',defaultVerbose,@isnumeric)
addParameter(p,'visualize',defaultVisualize,@isnumeric)
addParameter(p,'symmetryAngles',null(0), @isnumeric)

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
searchAngles = args.searchAngles;
filterAngles = symmetryAngles;
for sa = searchAngles(:)'
    filterAngles = [filterAngles mod(symmetryAngles+sa,2*pi)];
    filterAngles = unique(round(filterAngles,FLOAT_EQUALITY_PRECITION));
end
filterAngles = sort(squeeze(filterAngles));

%% Pregenerate a set of filtered images
if verbose > 0; disp('Generate set of filtered images');end;
if verbose > 0; tic; end;
filterSet = generateFilterSet(img, filterAngles, args.sigmas,visualize);
if verbose > 0; toc; end;

%% Scan angles
houghSet = cell(numel(symmetryAngles),1);
houghOffest = ceil(norm([imgWidth, imgHeight]));

if verbose > 0; disp('Start symmetry scanning'); end;
if verbose > 0; tic; end;

% If the Parallel computing toolbox is implemented, the parfor will run
% using a pool of worker threads. It should run like an ordinary for-loop
% otherwize.
parfor phiIndx = 1:1:numel(symmetryAngles)
    phi = symmetryAngles(phiIndx);
    houghAngle = phi - pi/2;
    
    % Finding the indices for the corresponding images in the set of pre
    % filtered images.
    if verbose > 0; fprintf('Phi: %.2fpi\n',phi/pi); end;     
    imageIndexes = zeros(2,size(searchAngles,2));
    for i = 1:size(searchAngles,2)
        theta = searchAngles(:,i);        
        theta1 = round(mod(phi+theta(1),2*pi),FLOAT_EQUALITY_PRECITION);
        theta2 = round(mod(phi+theta(2),2*pi),FLOAT_EQUALITY_PRECITION);        
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
    rotAngle= 180*(houghAngle)/pi;
    imgRot = imrotate(SIM,rotAngle,'bilinear');

    houghAngleOffset = min(0,round(imgHeight * sin(houghAngle)));        
    hough_line = zeros(1,2*houghOffest);      
    hough_line((1:size(imgRot,2)) + houghOffest + houghAngleOffset) = squeeze(sum(sum(imgRot,3),1));  
    houghSet{phiIndx} = hough_line;
        
end
if verbose > 0; toc; end;
%% Generate matrix for hough column

houghSpace = cat(1,houghSet{:});
[~, i] = max(houghSpace(:));
[phiIndx, rho]= ind2sub(size(houghSpace),i);

phi = symmetryAngles(phiIndx);
rho = rho - houghOffest;

% This version don't finds start and end positions.
lo = -1000;
hi = 1000;


%% Compute segments
segments = cell(1,1);
for i = 1:1
    theta = phi(i) - pi/2;
    [cx, cy] = pol2cart(theta,rho(i));
    [lx, ly] = pol2cart(theta+pi/2,lo(i));
    [hx, hy] = pol2cart(theta+pi/2,hi(i));

    segments{i} = [[cx+lx;cx+hx], [cy+ly;cy+hy]];

end


end






