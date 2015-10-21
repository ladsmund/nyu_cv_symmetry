function [rhos, values, lowerBounds, upperBounds] = selectCandidate(imgRot, parameters)

dWeight = ones(1,size(imgRot,3));

rhoHistogram = dWeight*squeeze(sum(imgRot,1))';

[peakValues, peakRhos] = findpeaks(rhoHistogram);
[~, sortIndices] = sort(peakValues(:),1,'descend');

indices = zeros(numel(sortIndices),1);
numberOfLines = 0;
for i = sortIndices'
    if sum(abs(peakRhos(indices(1:numberOfLines)) - peakRhos(i)) < parameters.miniumRhoDistance)
        continue
    end
    numberOfLines = numberOfLines + 1;
    indices(numberOfLines) = i;
    if numberOfLines >= parameters.numberOfLines
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
    
    energy = cumsum(lengthHistogram) / sum(lengthHistogram);
    lowerBounds(i) = find(energy > parameters.segmentHistogramQuantileLow,1);
    upperBounds(i) = find(energy > parameters.segmentHistogramQuantileHigh,1);   
end
end