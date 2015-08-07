function outputImg = outline2ellipse(img,embOutline,template)
% img: the original image
% embOutline: estimated embryo outline from the segmentation algorithm
% template: the ellipse template
% parts of the code were adapted from SPEX2
m = size(img,1);
n = size(img,2);

embRegion = poly2mask(embOutline(:,1),embOutline(:,2),m,n);


S = regionprops(embRegion, 'ConvexImage', 'BoundingBox');

embRegionConvex= S(1).ConvexImage;  % Convex-Hull smoothing
OutlineBoundingBox_one = max( ones(1,4), floor(S(1).BoundingBox));
% truncate to positive integers

S = regionprops(embRegionConvex, 'Orientation');
rotationAngle = -S(1).Orientation;
convexRegion = imrotate(embRegionConvex, rotationAngle, 'bilinear', 'loose');
imgTemp = img(OutlineBoundingBox_one(2) + (0:OutlineBoundingBox_one(4)-1), OutlineBoundingBox_one(1) + (0:OutlineBoundingBox_one(3)-1),:);  
imgTemp = imrotate(double(imgTemp), rotationAngle, 'bilinear', 'loose');

% now, find the min and max indices of the embryo region.
xHist = sum(convexRegion,1);
indTemp = find(xHist > 10);
xMin = min(indTemp); xMax = max(indTemp);
convexRegion = convexRegion(:,xMin:xMax);
imgTemp = imgTemp(:,xMin:xMax,:);
convexRegion = imresize(convexRegion, [size(convexRegion, 1), size(template,2)], 'nearest');
imgTemp = imresize(imgTemp, [size(imgTemp, 1), size(template,2)], 'nearest');% resize the region to match the width of the template


outputImg = zeros(size(template));
for j = 1 : size(template,2)
    if any(convexRegion(:,j))     
        outputImg(template(:,j)~=0,j,:) = imresize(imgTemp((convexRegion(:,j)~=0), j,:), [sum(template(:,j)), 1], 'nearest');
        % morphing the vertical dimension only...
    end;
end;
outputImg(~template) = 0;
