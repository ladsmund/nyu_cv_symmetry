function [SYM_rotated, houghRhoOffset, lengthOffset] = computeSymmetryMatrix(img, phi, parameters)
if nargin < 3
    parameters = get_default_parameters();
end
filterCombinator = parameters.filterCombinator;    
symmetryMetric = parameters.symmetryMetric;
searchAngles = parameters.searchAngles;
searchRange = parameters.searchRange;
sigmas = parameters.sigmas;
[imgHeight, imgWidth] = size(img);

biasFactor = parameters.distanceBiasAlpha * searchRange.^-2;
biasFactorMatrix = repmat(permute(biasFactor,[3,1,2]),imgHeight,imgWidth,1);

%% Compute Symmetry matrix
% The similarity matrix is a 3D matrix where the coordinates are the with,
% height and search distance.
SYM = zeros(imgHeight,imgWidth,length(searchRange));
for i = 1:numel(searchAngles)
    theta = searchAngles(i);    
    
    theta1 = mod(phi + theta - pi/2, 2*pi);
    theta2 = mod(phi+ pi/2 - theta,2*pi);

    J1 = morletFilter(img,theta1, sigmas);
    J2 = morletFilter(img,theta2, sigmas);

    for rIndx = 1:numel(searchRange)
        r = searchRange(rIndx);
        [dx, dy] = pol2cart(phi-pi/2,r);
        d = ([dx, dy]);
        J1t = imtranslate(J1,d);
        J2t = imtranslate(J2,-d);

        % Combining the two filtered images
        SYM(:,:,rIndx) = SYM(:,:,rIndx) + filterCombinator(J1t, J2t);
    end                
end

% Compute real values symmetry metric for voting    
SYM = symmetryMetric(SYM);
SYM = SYM + biasFactorMatrix;

%% Transform into rho/line/distance space
% By rotating the similarity matrix with respect to the houghAngle.
% The function imrotate used degree and rotates in the opposite direction.
    
houghAngle = phi - pi/2;
rotAngle= 180*(houghAngle)/pi;
SYM_rotated = imrotate(SYM,rotAngle,'bilinear');

houghRhoOffset = min(0,round(imgHeight * sin(houghAngle)));
lengthOffset = min(0,-round(imgWidth * sin(houghAngle)));
end