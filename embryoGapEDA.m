addpath('./registration')

% register embryo to template and run edge detector
imgPath = 'IM1/embImgs/debug_52.tif';
img = imread(imgPath);
registeredEmb = rgb2gray(registerEmbryo(img));
mask = registeredEmb==0; mask([1,end],:) = 1; mask(:,[1,end]) = 1;
mask = 1-imdilate(mask,strel('disk',50));
embEdges = edgeDetector(registeredEmb,1);


    
figure
imgScaled = imresize(embEdges,[300,600],'nearest');
maskScaled = imresize(mask,[300,600],'nearest');
%threshold edge detector responses
imgThresh(imgThresh<quantile(imgThresh(:),0.9)) = 0; 
imgFinal = imgThresh.*maskScaled;
imagesc(imgFinal)
