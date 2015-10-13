function [rhos, values, lowerBounds, upperBounds] = selectCandidate(imgRot, maxNumberOfLines, visualize)


if nargin < 3
    visualize = 0;
end

% dWeight = fspecial('gauss',[1,2*size(imgRot,3)],4);
% dWeight = dWeight(size(imgRot,3)+1:end);
dWeight = ones(1,size(imgRot,3));

% Generate histogram for rho values
% rhoHistogram = squeeze(sum(sum(imgRot,1),3));
rhoHistogram = dWeight*squeeze(sum(imgRot,1))';

rhoHistogramNorm = sqrt(rhoHistogram .* conj(rhoHistogram));
rhoHistogramNormBlur = conv(rhoHistogramNorm, fspecial('gaussian',[1,40],10),'same');

[peakValues, peakRhos] = findpeaks(rhoHistogramNormBlur);
[~, sortIndices] = sort(peakValues(:),1,'descend');

numberOfLines = min(maxNumberOfLines, numel(peakRhos));

rhos = peakRhos(sortIndices(1:numberOfLines));
values = peakValues(sortIndices(1:numberOfLines));

THRESSHOLD = 0.1;
lowerBounds = zeros(1,numberOfLines);
upperBounds = zeros(1,numberOfLines);
for i = 1:numberOfLines
    
    candidate = imgRot(:,rhos(i),:);

%     lengthHistogram = sum(candidate,3);
    lengthHistogram = squeeze(candidate) * dWeight';
    lengthHistogram = squeeze(sum(lengthHistogram,2));
    lengthHistogram_norm = sqrt(lengthHistogram .* conj(lengthHistogram));

    lengthHistogram_norm_mean = mean(lengthHistogram_norm);
    
    lowerBounds(i) = find((lengthHistogram_norm/lengthHistogram_norm_mean) > THRESSHOLD,1);
    upperBounds(i) = find((lengthHistogram_norm/lengthHistogram_norm_mean) > THRESSHOLD,1,'last');

end



if(visualize)   
    
    figure(200)
    plot(1:numel(rhoHistogram), rhoHistogramNorm)
    % plot3(1:numel(rhoHistogram), real(rhoHistogram), imag(rhoHistogram))
    hold on
    plot(1:numel(rhoHistogram), rhoHistogramNormBlur)
    % plot(1:numel(rhoHistogram)-1, 5000*zeroCross)
    hold off
    
    for i = 1:numberOfLines
        candidate = imgRot(:,rhos(i),:);

        cross_section = squeeze(sum(candidate,2));
        cross_section_norm = sqrt(cross_section .* conj(cross_section));

        subplot(1,numberOfLines,i);
        imagesc(cross_section_norm);
        title(rhos(i))
    end
    drawnow
%     waitforbuttonpress    
end



end