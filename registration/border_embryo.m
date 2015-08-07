function [bparam border img2t confidence] = border_embryo(inimg, interval, version, verbose)
% Segments the embryo
% inimg = RGB image
% interval = select 1:interval:end datapoints in result matrix
% version = version of finding border embryo (1-3)
% verbose = when set, say what's going on
% bparam = selected output parameters for finding embryo borders
% border = border coordinates of main object
% img2t = scaled image with variance
% confidence = confidence for assignments (0 is single embryo)
%   otherwise confidence is a bitfield
%       Bit     Meaning
%       1       embryos found and correctly separated
%       2       2 embryos touching the main one
%       3       Only one boundary found, second boundary is best guess
%       4       even number of boundaries
%       5       odd number of boundaries 

% parameters
Pvar_win = 3;
Pvar_thr = 2;
Pmo = 18;
Pbord_px = 50;
Pct_not_valid = 0.95;

if nargin < 2,
    verbose = 0;
end

confidence = 0;

cs = rgb2gray(inimg);
% imresize moved out for global consistency
% cs = imresize(cg, 0.5);

[nr,nc] = size(cs);


if verbose == 1, fprintf('Extracting the embryo with std_dev ...\n'); end;

% std deviation
B = im2col(cs,[Pvar_win,Pvar_win],'sliding');
B = double(B);
C = std(B);

% threshold it; all pixel with low variance to 0
Ct = C;
Ct(C < Pvar_thr) = 0;

% reconstitute image out of variance
img2t = zeros(nr, nc);
img2t(1:nr-Pvar_win+1,1:nc-Pvar_win+1) = reshape(Ct, nr-Pvar_win+1, nc-Pvar_win+1);


if verbose == 1, fprintf('Generating cleaned up binary image ...\n'); end;

% clean up image
im2l = logical(zeros(size(img2t)));
im2l(img2t >= 1) = 1;
es = bwmorph(im2l,'clean');
es = bwmorph(es,'dilate');
es = bwmorph(es,'majority',3);
bw = es;


if verbose == 1, fprintf('Generating outlines of center object ...\n'); end;

% take major object only
img3 = bwlabel(bw);
flagcenter = img3(round(nr/2-Pmo:nr/2+Pmo), round(nc/2-Pmo:nc/2+Pmo));
slist = unique(flagcenter(:));
for i=1:length(slist),
  lenlist(i) = nnz(flagcenter==slist(i));
end;
[slen, sidx] = sort(-lenlist);
rgnidx = slist(sidx(1));

% create bw image with major object
img3o = zeros(nr,nc);
img3o(img3 == rgnidx) = 1;
img3o = logical(img3o);

% find the boundaries
[B,L] = bwboundaries(img3o);
% Find the largest blob
boundary = []; A_max = 0;
A_notvalid = nc*nr*Pct_not_valid;

for i=1:length(B)
    A = false(nr, nc);
    A(sub2ind([nr, nc], B{i}(:,1), B{i}(:,2))) = 1;
    A_curr = regionprops(A,'FilledArea');
    if A_curr.FilledArea > A_notvalid,
        % seem to be some cases where the whole image is taken
        continue;
    end
    if A_curr.FilledArea > A_max,
        boundary = B{i};
        A_max = A_curr.FilledArea;
    end
end

x = boundary(:, 2) ./ nc;
y = boundary(:, 1) ./ nr;
nb = length(boundary);

% calculate angle & distance (currently from image center)
yd = y - ((nr/2) / nr);
xd = x - ((nc/2) / nc);
d = sqrt(xd.^2 + yd.^2);
a = atan2(yd, xd);

% find variance in window
%da_cont = [d, a, boundary(:,2), boundary(:,1)];
da_cont = [boundary(:,2), boundary(:,1), x, y, a, d ];
da_cont(:,7) = 1: length(da_cont);
da_conts = sortrows(da_cont,5);
dev = im2col([da_conts(:,6); da_conts(1:4,6)],[5 1], 'sliding');
dev_std = std(dev);
da_conts(:,8) = dev_std;
% smoothen variance
%dev_std_csm = im2col(dev_std,[1,10],'sliding');
%dev_std_sm = mean(dev_std_csm);

% 12/1 use non-normalized version
border = da_conts;
% border = [x y];


if verbose == 1, fprintf('Testing for multiple embryos ... '); end;

% find if there are multiple embryos
bord_px = zeros(4,1);
bord_px(1) = length(find(border(:,2) == 1));
bord_px(2) = length(find(border(:,2) > nr - 1));
bord_px(3) = length(find(border(:,1) == 1));
bord_px(4) = length(find(border(:,1) > nc - 1));

membryo  = 0;
if max(bord_px) >= Pbord_px
    if verbose == 1, fprintf('found\n'); end;
    membryo = 1;
else
    if verbose == 1, fprintf('not found\n'); end;
end


switch(version)
    case 1
        % version 1
        % interesting parameters are 6:7
        % 6: distance from center
        % 7: std-dev in region
        bparam = [da_conts(1:interval:nb,:) (dev_std(1:interval:nb))'];

    case 2
        % version 2
        % interesting parameters are 7:8
        % 7: std-dev in region
        % 8: minimum distance of pixel from boundary
        da_dist(1,:) = da_cont(:,3);
        da_dist(2,:) = abs(da_cont(:,3) - 760);
        da_dist(3,:) = da_cont(:,4);
        da_dist(4,:) = abs(da_cont(:,4) - 540);
        da_dist(5,:) = 1:length(da_dist);

        da_conts = sortrows(da_conts, 7);
        bparam = [ da_conts(1:interval:nb,1:6) da_conts(1:interval:nb,8) min(da_dist(:, 1:interval:nb))' ];
    
    case 3
        % version 3
        if verbose == 1, fprintf('Finding overlaps for border lines (version 3) ...\n'); end;
        da_conts(:,9) = fix(da_conts(:,5) .* (180/pi));
        coord_lines = [];
        for i = min(da_conts(:,9)):max(da_conts(:,9)),
            coord = da_conts(find(da_conts(:,9) == i),1:2);
            dev = da_conts(find(da_conts(:,9) == i),8);
            dev = median(dev);
            nlines = 0;
            if size(coord,1) > 1
                coord = coord - repmat((min(coord) - 1), length(coord),1);
                coord_img = zeros(max(coord));
                coord_img(sub2ind(size(coord_img), coord(:,1), coord(:,2))) = 1;
                if (length(unique(coord_img))) == 1
                    nlines = 1;
                else
                    nlines = length(unique(bwlabel(coord_img)))-1;
                end
            else
                nlines = 1;
            end
            coord_lines = [coord_lines; i nlines dev];
        end
        bparam = coord_lines;
        
        if membryo == 0
            % if no multiple embryos are found, skip the rest
            bparam = [ bparam zeros(length(coord_lines),1) ];
            return;
        end
                
        if verbose == 1, fprintf('Pinpointing the embryo overlaps ...\n'); end;
        % calculate threshold for std_dev of distances for touching embryos
        coord_sdev = sort(coord_lines(:,3));
        touch_dmax = max(abs(diff(coord_sdev))) * 0.5;
        % touch_dsort = sort(abs(diff(coord_sdev)));
        % touch_dmax = touch_dsort(length(coord_sdev) - 5);
        touch_dthresh = find(coord_sdev >= touch_dmax);
        touch_max = coord_sdev(touch_dthresh(1));
        if verbose == 1, touch_max, end;
        %touch_min = coord_sdev(200);
        
        % assign everything with more than 3 non-touching lines to a block
        touch = find(coord_lines(:,2) >= 2); % 11/29 test with 2 to avoid dicontinuities
        touch_c = find(diff(touch) > 1);
        touch_block = [touch(1) 0 0];
        for i = 1:length(touch_c)
            touch_block(i, 2) = touch(touch_c(i));
            touch_block = [ touch_block; touch(touch_c(i)+1) 0 0];
        end
        touch_block(size(touch_block,1), 2) = touch(end);
        
        % calculate the std_dev for each block and find those below
        % threshold to be deleted
        touch_del = [];
        for i = 1:size(touch_block,1)
            touch_block(i,3) = max(coord_lines(touch_block(i,1):touch_block(i,2),3));
            % touch_block(i,3) = mean(coord_lines(touch_block(i,1):touch_block(i,2),3));
            if touch_block(i,3) < touch_max,
                touch_del = [ touch_del i ];
            end
        end
        if verbose == 1, touch_block, end;
        
        % delete entries below threshold and (11/29) more than 2 entries found
        if size(touch_block,1) > 2
            % save sorted touch_block 
            touch_sblock = sortrows(touch_block,3);
            if length(touch_del) > 0
                % delete rest
                touch_block(touch_del,:) = [];
            end
        end
        if verbose == 1, touch_block, end;
        
        confidence = 0;
        em = zeros(length(coord_lines),1);        
        touch_no = size(touch_block,1);
        if touch_no < 2,
            em = assign_em(touch_sblock(end-1:end,:), em);
            if verbose == 1, fprintf('Uncertain (less than 2 touchpoints: %d)!\n', touch_no); end;
            confidence = bitset(confidence, 3);
            confidence = bitset(confidence, 5);
        elseif mod(touch_no, 2) == 0
            if touch_no > 4
                em = assign_em(touch_sblock(end-1:end,:), em);
                if verbose == 1, fprintf('Uncertain (more than 2 embryos: %d) - taking one best guess embryo!\n', touch_no); end;
                confidence = bitset(confidence, 4);
            else
                em = assign_em(touch_block, em);
                confidence = bitset(confidence, 1);
                if touch_no == 4
                    confidence = bitset(confidence, 2);
                end
            end
        else
            em = assign_em(touch_sblock(end-1:end,:), em);
            if verbose == 1, fprintf('Uncertain (odd number %d) - taking one best guess embryo!\n', touch_no); end;
            confidence = bitset(confidence, 5);
        end
        bparam = [bparam em];
end


function [em] = assign_em(touch_block, em_in)
% assigns "emmisions" to angles (1 = overlapping embryo)

em = em_in;
touch_block = sortrows(touch_block, 1);

touch_no = size(touch_block,1);
touch_block(:,3) = circshift(touch_block(:,2), 1);
touch_block(:,4) = touch_block(:,1) - touch_block(:,3);
touch_block(1,4) = touch_block(1,4) + 360;
touch_block = sortrows(touch_block,4);
touch_block = touch_block(1:touch_no / 2,:);
for i = 1:size(touch_block,1)
    if touch_block(i,3) < touch_block(i,1)
        em(touch_block(i,3):touch_block(i,1)) = i;
    else
        em(1:touch_block(i,1)) = i;
        em(touch_block(i,3):length(em)) = i;
    end
end
