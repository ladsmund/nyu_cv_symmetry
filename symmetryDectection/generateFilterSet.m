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
        filterSet.filtered{i,j} = imgPadFilterd(xPad+1:end-xPad, yPad+1:end-yPad);
        
%         if visualize >= 1
%             fig = figure(42);
%             
%             subplot(2, length(filterAngles),j)
%             imagesc(imag(filterSet.filtered{i,j}))
%             [px, py] = pol2cart(filterAngles(j),20);
%             line(xPad/2+[0, px], yPad/2 + [0, py])
%             title(filterAngles(j))
%             
%             subplot(2,length(filterAngles),j+length(filterAngles))
%             imagesc(real(filterSet.filtered{i,j}))
%             [px, py] = pol2cart(filterAngles(j),20);
%             line(xPad/2+[0, px], yPad/2 + [0, py])
%             title(filterAngles(j))
%             
            
%             drawnow
%             subplot(2,2,1)
%             imagesc(imag(kernel))
%             [px, py] = pol2cart(filterAngles(j),20);
%             line(xPad/2+[0, px], yPad/2 + [0, py])
%             title(filterAngles(j))
%             subplot(2,2,2)
%             imagesc(imag(filterSet.filtered{i,j}))
%             [px, py] = pol2cart(filterAngles(j),20);
%             line(xPad/2+[0, px], yPad/2 + [0, py])
%             title(filterAngles(j))
%             subplot(2,2,3)
%             imagesc(real(kernel))
%             line(xPad/2+[0, px], yPad/2 + [0, py])
%             title(filterAngles(j))
%             subplot(2,2,4)
%             imagesc(real(filterSet.filtered{i,j}))
%             [px, py] = pol2cart(filterAngles(j),20);
%             line(xPad/2+[0, px], yPad/2 + [0, py])
%             title(filterAngles(j))
%             drawnow
%         end
        
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

