function filtered = morletFilter(img, angle, sigma)
% Create kernel
kernel = morletWavelet(sigma,angle);  
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

end