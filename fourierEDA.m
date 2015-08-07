%loadPath = './IM1/embImgs/';
loadPath = './IM5/javaData/';
filesTemp = dir([loadPath,'debug*']); 
nFiles = length(filesTemp);

width = 840; height = 631;
%width = 1232; height = 816;
r = sqrt((height/2)^2+(width/2)^2);
propRadius = 0.10:0.0025:0.17;
nbhdBuffer = 0.02;
nRadii = length(propRadius);
densities = zeros(nRadii-4,nFiles);

fileNames = cell(nFiles,1);
imgScale = zeros(nRadii-4,nFiles);
for j=1:nFiles
    
	fprintf('img: %d\n', j)
    fileNames{j} = filesTemp(j).name;
    %img = imread([loadPath,filesTemp(j).name]);
    %fftImg = fftshift(fft2(img));
    fftImg = importdata([loadPath,filesTemp(j).name]);
    if length(fftImg)~=530040
        fprintf('img: %d incorrect size\n', j)
        continue
    end
    fftImg = reshape(fftImg, 631, 840);
    rootName = strsplit(filesTemp(j).name,'.'); rootName = rootName(1);

    for i = 1:(nRadii-4)
	
        nbhdMax = propRadius(i);
        nbhdMin = propRadius(i)-nbhdBuffer;
        bpfRadiusMax = propRadius(i+4)*r; bpfRadiusMin = propRadius(i)*r;
        nbhdRadiusMax = nbhdMax*r; nbhdRadiusMin = nbhdMin*r;
        bpfArea = pi*bpfRadiusMax^2-pi*bpfRadiusMin^2;
        bpf = zeros(size(fftImg));

        baseY = linspace(-height/2,height/2,height);
        baseX = linspace(-width/2,width/2,width);
        [y,x] = meshgrid(baseX, baseY);
        bpf(x.^2+y.^2<bpfRadiusMax^2) = 1;
        bpf(x.^2+y.^2<bpfRadiusMin^2) = 0;

        bpfNbhd = zeros(size(fftImg));
        bpfNbhd(x.^2+y.^2<bpfRadiusMax^2) = 1;
        bpfNbhd(x.^2+y.^2<nbhdRadiusMin^2) = 0;
        bpfNbhdArea = pi*bpfRadiusMax^2-pi*nbhdRadiusMin^2;

        % cut out some proportion of low frequencies
        midCut = 0.0075; 
        xThresh = height*midCut;
        yThresh = width*midCut;
        xLow = -x<xThresh; xHigh = x<xThresh;
        yLow = -y<yThresh; yHigh = y<yThresh;
			   
        bpf(intersect(find(xLow),find(xHigh))) = 0;
        bpf(intersect(find(yLow),find(yHigh))) = 0;
        
        filteredImg = log(abs(fftImg)).*bpf;
        densitiesTemp = sum(filteredImg(:))/bpfArea;
        filteredImgNbhd = log(abs(fftImg)).*bpfNbhd;
        imgScaleTemp = sum(filteredImgNbhd(:))/bpfNbhdArea;
        densities(i,j) = densitiesTemp/imgScaleTemp;
    end
    
    
    %{
    slopes = linspace(-1,1,20);
    nSlopes = length(slopes);
    features = cell(nSlopes,1);
    for k=1:nSlopes
        
        buff = 5;
        xIdx = -(width/2):1:(width/2);
        yIdx = floor(slopes(k).* xIdx);
        xIdx = xIdx+(width/2);
        yIdx = yIdx+(height/2);
        yLow = yIdx - buff; yHigh = yIdx + buff;
        nLin = length(xIdx);
        tempFeatures = zeros(nLin,1);
        for i=1:nLin
            yLowTemp = max(yLow(i),1);
            yHighTemp = min(yHigh(i),height);
            yRan = yLowTemp:yHighTemp;
            nRan = length(yRan);
            xRan = xIdx(i)*ones(nRan,1)';
            if xIdx(i)~=0
                linIdcs = sub2ind(size(filteredImg),yRan,xRan);
                tempFeatures(i) = mean(filteredImg(linIdcs));
            else
                tempFeatures=NaN;
            end
        end
        features{k} = tempFeatures;
    end
    features = cell2mat(features);
    
    
    
    for i = 1:nSlopes
      minRad = cos(atan(slopes(i)))*bpfRadiusMin;
      maxRad = cos(atan(slopes(i)))*bpfRadiusMax;
      f = figure('Visible','off');
      plot(features(i,:)); 
      hold on
      plot(width/2+minRad*ones(1,11),0:10)
      plot(width/2+maxRad*ones(1,11),0:10)
      plot(width/2-minRad*ones(1,11),0:10)
      plot(width/2-maxRad*ones(1,11),0:10)
      print(f,[imgDir,cell2mat(rootName),'_',int2str(i)],'-dpng');
    end
    %}


%{
    f = figure('Visible','off')
  				       
    yGridL = y(abs(y)<bpfRadiusMin); yGridH = y(abs(y)<bpfRadiusMax);
    xC = width/2; yC = height/2;
    plot(yGridL+xC,sqrt(bpfRadiusMin^2-yGridL.^2)+yC); plot(yGridL+xC,-sqrt(bpfRadiusMin^2-yGridL.^2)+yC);
    plot(yGridH+xC,sqrt(bpfRadiusMax^2-yGridH.^2)+yC); plot(yGridH+xC,-sqrt(bpfRadiusMax^2-yGridH.^2)+yC);
    print(f,[imgDir,cell2mat(rootName)],'-dpng')
      
%}
end

%{
f = figure('Visible','off')
for j=20:nFiles
plot(densities(:,j));
	hold on
end
	legend(int2str(20:nFiles))
	print(f,[imgDir,'plot'],'-dpng')
	  
%}
plot(densities,':b')

i = find(ismember(fileNames,'debug_194.txt'));
hold on
plot(densities(:,i),'g')

score = densities;
[val,idx] = sort(score);
fprintf('Opt Img: %s', fileNames{idx(end)})



%{
imgIdcs = regexp(fileNames,'[_.]','split');
imgIdcs = cellfun(@(c)c(2),imgIdcs);
imgIdcs = cellfun(@str2num,imgIdcs);
minImgIdx = min(imgIdcs);
maxImgIdx = max(imgIdcs);


f = figure('Visible','off');
[vals,order] = sort(imgIdcs);
plot(vals,densities(1,order)./imgScale(order))
xlabel('Img index')
ylabel('score')

[s,i] = sort(densities(1,:)./imgScale );%+ (densities(1,:)./imgScale2)*0.5);
fprintf('Max Score: %s ',fileNames{i(end)})    

hold on
maxIdx = (minImgIdx-1)+i(end);
maxVal = s(end);
text(maxIdx,maxVal,int2str(maxIdx));

maxIdx2 = (minImgIdx-1)+i(end-1);
maxVal2 = s(end-1);
text(maxIdx2,maxVal2,int2str(maxIdx2));

print(f,[imgDir,'scores'],'-dpng')
%}
