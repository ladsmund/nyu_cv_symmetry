function show_result(data, image_path)
%%

if nargin < 2
    image_path = data.image_path;
end

img = imread(image_path);



%%

imshow(img)
colors = hot(numel(data.segments));
hold on

segment = data.segments{1};
colorIndx = numel(data.segments);
for i = 1:numel(data.segments)
    segment = data.segments{i};
    line([segment(1,1),segment(2,1)], [segment(1,2),segment(2,2)],'Color',colors(colorIndx,:), 'LineWidth', 2)        
    colorIndx = colorIndx - 1;
end

segment = data.segments{1};
line([segment(1,1),segment(2,1)], [segment(1,2),segment(2,2)],'Color',[1,1,1], 'LineWidth', 4)        
line([segment(1,1),segment(2,1)], [segment(1,2),segment(2,2)],'Color',[0,0,0], 'LineWidth', 4, 'LineStyle','--')


hold off

% frame = getframe(fig);
% imwrite(frame.cdata,[outputFolder, name, ext])

end