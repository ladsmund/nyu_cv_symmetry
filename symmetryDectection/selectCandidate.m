function [rhos, values, lowerBounds, upperBounds] = selectCandidate(imgRot, maxNumberOfLines, sigma, visualize)


if nargin < 4
    visualize = 0;
end

dWeight = ones(1,size(imgRot,3));

rhoHistogram = dWeight*squeeze(sum(imgRot,1))';

rhoHistogramNorm = sqrt(rhoHistogram .* conj(rhoHistogram));

[peakValues, peakRhos] = findpeaks(rhoHistogramNorm);
[~, sortIndices] = sort(peakValues(:),1,'descend');

indices = zeros(numel(sortIndices),1);
numberOfLines = 0;
for i = sortIndices'
    if sum(abs(peakRhos(indices(1:numberOfLines)) - peakRhos(i)) < 3*sigma)
        continue
    end
    numberOfLines = numberOfLines + 1;
    indices(numberOfLines) = i;
    if numberOfLines >= maxNumberOfLines
        break;
    end
end
sortIndices = indices(1:numberOfLines);

rhos = peakRhos(sortIndices(1:numberOfLines));
values = peakValues(sortIndices(1:numberOfLines));

lowerBounds = zeros(1,numberOfLines);
upperBounds = zeros(1,numberOfLines);
for i = 1:numberOfLines
    
    candidate = imgRot(:,rhos(i),:);
    lengthHistogram = squeeze(candidate) * dWeight';
    lengthHistogram = squeeze(sum(lengthHistogram,2));
    lengthHistogram_norm = sqrt(lengthHistogram .* conj(lengthHistogram));
    
    energy = cumsum(lengthHistogram_norm) / sum(lengthHistogram_norm);
    lowerBounds(i) = find(energy > .02,1);
    upperBounds(i) = find(energy > .98,1);
    
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