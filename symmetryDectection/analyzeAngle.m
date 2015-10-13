function [lines] = analyzeAngle(symAngle, imageIndexes, filterSet, searchRange, maxNumberOfLines)

houghAngle = symAngle - pi/2;

imgHeight= size(filterSet.img,1);
imgWidth = size(filterSet.img,2);

% The similarity matrix is a 3D matrix where the coordinates are the with,
% height and search distance.
SYM = zeros(imgHeight,imgWidth,length(searchRange));

%% Compute similarity matrix
for i = imageIndexes
    tIndx1 = i(1);
    tIndx2 = i(2);

    J1 = filterSet.filtered{tIndx1};
    J2 = filterSet.filtered{tIndx2};
    
    for rIndx = 1:numel(searchRange)
        r = searchRange(rIndx);
        [dx, dy] = pol2cart(symAngle-pi/2,r);
        d = ([dx, dy]);        
        J1t = imtranslate(J1,d);
        J2t = imtranslate(J2,-d);
        
%         symmetryShow(J1t,J2t,r);
        % TODO: SYM might has to be squared because it's a sum om squared
        % values
%         SYM(:,:,rIndx) = SYM(:,:,rIndx) + sqrt(real(J1t .* conj(J2t)));
%         SYM(:,:,rIndx) = SYM(:,:,rIndx) + (real(J1t .* conj(J2t)));
%         SYM(:,:,rIndx) = SYM(:,:,rIndx) + sqrt(J1t .* conj(J2t));
        SYM(:,:,rIndx) = SYM(:,:,rIndx) + (J1t .* conj(J2t));
%         SYM(:,:,rIndx) = SYM(:,:,rIndx) + sqrt(J1t + J2t);

    end                
end


% SYM = sqrt(real(SYM .* conj(SYM)));
SYM = real(sqrt(SYM .* conj(SYM)));
% SYM = real(SYM);

% TODO: Add a gaussian bias related to distance because people usualy
% prefere small distances when looking for symmetry.

%% Transform into rho/line/distance length space
% By rotating the similarity matrix with respect to phi.

% The function imrotate used degree and rotates in the opposite direction.
rotAngle= 180*(houghAngle)/pi;
imgRot = imrotate(SYM,rotAngle,'bilinear');

houghRhoOffset = min(0,round(imgHeight * sin(houghAngle)));

lenOffset = min(0,-round(imgWidth * sin(houghAngle)));
% lenOffset = 0;
% lengthOffset = min(0,round(imgHeight * cos(houghAngle)));
lengthOffset = 0;

[rhos, values, lowerBounds, upperBounds] = selectCandidate(imgRot,maxNumberOfLines,0);

rhos = rhos + houghRhoOffset;
lowerBounds = lowerBounds + lengthOffset + lenOffset;
upperBounds = upperBounds + lengthOffset + lenOffset;

lines = [rhos; values; lowerBounds; upperBounds];
end