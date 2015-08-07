function registered = registerEmbryo(img)

[m,~,~] = size(img);
fac = 900/m;

% rescale the image
img = imresize(img,fac);
img = cat(3,img,img,img);
[height,width] = size(img);
template = generateTemplate(width,height,3);

% find the boundary of the embryo
[x,y,~] = fembryo(img,3,0);

% register the embryo on an ellipse
registered = uint8(outline2ellipse(img,[x,y],template));
