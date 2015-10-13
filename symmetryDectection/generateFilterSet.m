function filterSet = generateFilterSet(img, thetas, sigmas, visualize)

if nargin < 4
    visualize = 0;
end

filterSet = struct('img',img,'thetas',thetas,'sigmas',sigmas);
filterSet.filtered = cell(length(sigmas), length(thetas));

for i = 1:length(sigmas)
    for j = 1:length(thetas)
        kernel = morletWavelet(sigmas(i),thetas(j));
        
        xPad = size(kernel,1);
        yPad = size(kernel,2);
        
        % Padd
        imgPad = padarray(img,[xPad yPad],'replicate');
        imgPadFilterd = conv2(imgPad, kernel, 'same');
        filterSet.filtered{i,j} = imgPadFilterd(xPad+1:end-xPad, yPad+1:end-yPad);
        
        if visualize >= 1
            fig = figure(42);
            subplot(2,2,1)
            imagesc(imag(kernel))
            [px, py] = pol2cart(thetas(j),20);
            line(xPad/2+[0, px], yPad/2 + [0, py])
            title(thetas(j))
            subplot(2,2,2)
            imagesc(imag(filterSet.filtered{i,j}))
            [px, py] = pol2cart(thetas(j),20);
            line(xPad/2+[0, px], yPad/2 + [0, py])
            title(thetas(j))
            subplot(2,2,3)
            imagesc(real(kernel))
            line(xPad/2+[0, px], yPad/2 + [0, py])
            title(thetas(j))
            subplot(2,2,4)
            imagesc(real(filterSet.filtered{i,j}))
            [px, py] = pol2cart(thetas(j),20);
            line(xPad/2+[0, px], yPad/2 + [0, py])
            title(thetas(j))
            drawnow
        end
        
    end
end
if visualize >= 2
    close(fig);
end
end