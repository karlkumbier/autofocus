function convScore = edgeDetector(img, nAngles)

filtWidth = 15;
filtSig = 2;

angles = linspace(0,2*pi,nAngles+1);
angles = angles(1:(end-1));

% The spatial component (the same at every scale)
filtX = fspecial('gaussian', [1 filtWidth], filtSig);
filtY = fspecial('gaussian', [1 filtWidth], filtSig)';

% Filter image with the 2d filter
movF = convn(convn(img, filtX, 'same'), filtY, 'same');

% Now take the gradient
[I_X, I_Y] = gradient(movF);

% Loop through angles and get the directional gradients
alphas = zeros(size(I_X,1),size(I_X,2), nAngles);
meanResp = zeros(nAngles,1);
for j = 1:nAngles

    theta = angles(j);
    alphaFull =  max(cos(theta) * I_X + sin(theta) * I_Y,0);
    alphas(:,:,j) = alphaFull;
    meanResp(j) = mean2(alphaFull);
end

convScore = log(sum(alphas(:)));

