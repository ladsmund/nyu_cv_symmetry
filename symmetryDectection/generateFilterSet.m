function filterSet = generateFilterSet(img, filterAngles, sigmas, visualize)

if nargin < 4
    visualize = 0;
end

filterSet = struct('img',img,'thetas',filterAngles,'sigmas',sigmas);
filterSet.filtered = cell(length(sigmas), length(filterAngles));

for i = 1:length(sigmas)
    for j = 1:length(filterAngles)
        kernel = morletWavelet(sigmas(i),filterAngles(j));
        
        xPad = size(kernel,1);
        yPad = size(kernel,2);
        
        % Padd
        imgPad = padarray(img,[xPad yPad],'replicate');
        imgPadFilterd = conv2(imgPad, kernel, 'same');
        filtered = imgPadFilterd(xPad+1:end-xPad, yPad+1:end-yPad);

        % Use only central square of image
        [v, minIndx] = min(size(filtered));
        d = abs(diff(size(filtered)));
        if minIndx == 1
            filtered(:,1:floor(d/2)) = 0;
            filtered(:,(end-ceil(d/2)):end) = 0;
        else
            filtered(1:floor(d/2),:) = 0;
            filtered((end-floor(d/2)):end,:) = 0;
        end
        
        filterSet.filtered{i,j} = filtered;
        
        
    end
end

if visualize >= 1
    visualize_filter_set(filterSet);
end
if visualize >= 2
    close(fig);
end
end

function visualize_filter_set(filterSet)

fig = figure(42);
xPad = 50;
yPad = 50;
number_of_filters = length(filterSet.thetas);

max_real = 0;
max_imag = 0;
for j = 1:number_of_filters
    max_real = max([max_real, max(max(abs(real(filterSet.filtered{1,j}))))]);
    max_imag = max([max_imag, max(max(abs(imag(filterSet.filtered{1,j}))))]);
end

for j = 1:number_of_filters

    subplot(2, number_of_filters,j)
    imshow(imag(filterSet.filtered{1,j}) / max_imag / 2 + .5)
    [px, py] = pol2cart(filterSet.thetas(j),20);
    line(xPad/2+[0, px], yPad/2 + [0, py])
    title(filterSet.thetas(j))

    subplot(2,number_of_filters,j+number_of_filters)
    imshow(real(filterSet.filtered{1,j} / max_real / 2 + .5))
    [px, py] = pol2cart(filterSet.thetas(j),20);
    line(xPad/2+[0, px], yPad/2 + [0, py])
    title(filterSet.thetas(j))
    
end

drawnow

end

