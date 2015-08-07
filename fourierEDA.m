nDir = 24;
optImgs = cell(nDir,1);
width = 1232; height = 816;
radius = sqrt((height/2)^2+(width/2)^2);
filterSizes = [0.1,0.125];
nbhdBuffer = 0.02;
nFilters = length(filterSizes);

for k = 1:nDir
    
    fprintf('processing directory: %d\n', k)
    loadPath = ['./IM',int2str(k),'/embImgs/'];
    filesTemp = dir([loadPath,'debug*']);
    nFiles = length(filesTemp);
    
    
    filterScores = zeros(nFilters-1,nFiles);
    % calculate magnitude of fourier transform in regions defined by
    % filterSizes
    fileNames = cell(nFiles,1);
    for j=1:nFiles
        
        fileNames{j} = filesTemp(j).name;
        img = imread([loadPath,filesTemp(j).name]);
        if size(img,1) ~= height
            continue
        end
        fftImg = fftshift(fft2(img));
        
        for i = 1:(nFilters-1)
            
            
            bpfRadiusMax = filterSizes(i+1)*radius;
            bpfRadiusMin = filterSizes(i)*radius;
            nbhdRadiusMax = filterSizes(i)*radius;
            nbhdRadiusMin = (filterSizes(i)-nbhdBuffer)*radius;
            bpfArea = pi*bpfRadiusMax^2-pi*bpfRadiusMin^2;
            bpfNbhdArea = pi*bpfRadiusMax^2-pi*nbhdRadiusMin^2;
            
            baseY = linspace(-height/2,height/2,height);
            baseX = linspace(-width/2,width/2,width);
            [y,x] = meshgrid(baseX, baseY);
            
            % set bpf for region
            bpf = zeros(size(fftImg));
            bpf(x.^2+y.^2<bpfRadiusMax^2) = 1;
            bpf(x.^2+y.^2<bpfRadiusMin^2) = 0;
            
            % set bpf for local reweighting
            bpfNbhd = zeros(size(fftImg));
            bpfNbhd(x.^2+y.^2<bpfRadiusMax^2) = 1;
            bpfNbhd(x.^2+y.^2<nbhdRadiusMin^2) = 0;
            
            
            % cut out some proportion of low frequencies
            midCut = 0.0075;
            xThresh = height*midCut;
            yThresh = width*midCut;
            xLow = -x<xThresh; xHigh = x<xThresh;
            yLow = -y<yThresh; yHigh = y<yThresh;
            
            bpf(intersect(find(xLow),find(xHigh))) = 0;
            bpf(intersect(find(yLow),find(yHigh))) = 0;
            
            filteredImg = log(abs(fftImg)).*bpf;
            scoresTemp = sum(filteredImg(:))/bpfArea;
            filteredImgNbhd = log(abs(fftImg)).*bpfNbhd;
            imgRescaleTemp = sum(filteredImgNbhd(:))/bpfNbhdArea;
            filterScores(i,j) = scoresTemp/imgRescaleTemp;
        end
        
        
    end
    


    [val,idx] = sort(filterScores);
    fprintf('Opt Img: %s \n', fileNames{idx(end)})
    optImgs{k} = [loadPath,': ',fileNames{idx(end)}];

end


