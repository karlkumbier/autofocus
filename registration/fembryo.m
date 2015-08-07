function [x y score] = fembryo(inimg, version, verbose)
% Segments the embryo

Param_SNalpha = 0.5;
Param_SNgamma = 2;
Param_SNiter = 3;
Pbord_px = 30;

bverbose = 0;
if verbose > 0, 
    fprintf('Finding embryo border ...\n');
    bverbose = 1;
    if verbose > 1, figure; end;
end;

s_touch = 0;

[bparam, border, img2t, score] = border_embryo(inimg,1,version,bverbose);

[nr,nc] = size(img2t);

% vcanvas = zeros(nr,nc);
% vcanvas(sub2ind([nr, nc], border(:,2), border(:,1))) = 1;

tstate = bparam(:,4);

b_angle = fix(border(:,5) .* (180/pi)) + 180;

if length(find(tstate ~= 0)) > 0
    % touching embryos detected
    if verbose > 0, fprintf('Detected multiple embryos, interpolating main embryo ...\n'); end;
    s_touch = 1;
    score = bitshift(score, 2);

    if verbose > 1
        % start plot for splines with the border
        subplot(3,1,1);
        hold on; 
        plot(border(:,1), border(:,2)); 
    end
    
    % storage for calculated spline coordinates
    sp_xsp = [];
    sp_ysp = [];
    sp_range = zeros(360,1);
    
    for touch_no = 1:max(unique(tstate))

        touch = find(tstate == touch_no);

        touch_c(1) = touch(1);
        touch_c(2) = touch(end);

        if touch(1) == 1 && touch(end) == 360
            c = find(diff(touch) > 1);
            touch_c(2) = touch(c(1));
            touch_c(1) = touch(c(1)+1);
        end

        sp_t = [ touch_c(1) - 10, touch_c(1) - 2, touch_c(2) + 2, touch_c(2) + 10 ];
        if verbose > 0, sp_t, end;

        c = find(sp_t < 0);
        if length(c) >= 0 
            sp_t(c) = sp_t(c) + 360;
        end
           
        max_sp_loop = 5;
        sp_loop = 1;
        while 1
            % until we find a good spline we'll loop
            sp_x = [];
            sp_y = [];
            for i=1:length(sp_t)
                m_angle = find(b_angle == sp_t(i));
                % border(m_angle,:)
                [m_min, m_pos] = min(border(m_angle, 6)); %TODO: sanity check if min is in line with close values
                m_pos = m_pos + min(m_angle) - 1;
                x = border(m_pos,1);
                y = border(m_pos,2);
                sp_x = [sp_x, fix(x)];
                sp_y = [sp_y, fix(y)];
            end

            % for spline generation make sure that values are monotonous
            if all(diff(sp_y) >= 0) || all(diff(sp_y) <= 0)
                sp_mon = sp_y;
            else
                sp_mon = sp_x;
            end
            if any(diff(sp_mon) == 0)
                % if they are not monotonous, subtract/add 3 degrees from angle
                % and try again
                sp_loop = sp_loop + 1;
                if sp_loop > max_sp_loop
                    error('No suitable spline found');
                end
                if verbose > 0, fprintf('Not monotous, adjusting angles\n'); end;
                sp_monpos = find(diff(sp_mon) == 0);
                if sp_monpos(1) < 3
                    sp_t = sp_t + [-3 -3 0 0];
                else
                    sp_t = sp_t + [0 0 3 3];
                end
            else
                break;
            end     
        end

        sp_xy = [sp_x; sp_y]';

        % make the spline either in x or y (dependent on which one is
        % monotonous)
        if verbose > 0, fprintf('Making spline: '); end;
        if all(diff(sp_y) > 0) || all(diff(sp_y) < 0)
            if verbose > 0, fprintf('Monoticity in y -> spline along y axis\n'); end;
            sp_xy = sortrows(sp_xy, 2); sp_x = sp_xy(:,1)'; sp_y = sp_xy(:,2)';
            sp_s = [min(sp_y):max(sp_y)];
            sp_c = spline(sp_y, sp_x, sp_s);
            sp_ysp = [sp_ysp, sp_s];
            sp_xsp = [sp_xsp, sp_c];
        elseif all(diff(sp_x) > 0) || all(diff(sp_x) < 0)
            if verbose > 0, fprintf('Monoticity in x -> spline along x axis\n'); end;
            sp_xy = sortrows(sp_xy, 1); sp_x = sp_xy(:,1)'; sp_y = sp_xy(:,2)';
            [sp_x; sp_y]
            sp_s = [min(sp_x):max(sp_x)];
            sp_c = spline(sp_x, sp_y, sp_s);
            sp_xsp = [sp_xsp sp_s];
            sp_ysp = [sp_ysp sp_c];
        else
            if verbose > 0, fprintf('No monoticity in x or y found - error creating spline\n'); end;
            continue;
        end

        if verbose > 1
            % plot the spline
            plot(sp_x,sp_y,'o',sp_xsp,sp_ysp)
        end

        % more complicated version than just find(sp_angle == i) in case there
        % are missing data in sp_angle
        if sp_t(2) > sp_t(3)
            sp_range(sp_t(2):360) = 1;
            sp_range(1:sp_t(3)) = 1;
        else
            sp_range(sp_t(2):sp_t(3)) = 1;
        end
    end
 
    if verbose > 1
        % finish the plot with the splines
        hold off;
    end
    
    % calculate angles for spline
    sp_angle = fix(atan2((sp_ysp ./ nr) - (1/2), (sp_xsp ./ nc) - (1/2)) .* (180/pi)) + 180;

    % take point for every angle excluding the +-2 arc in sp_t as outline
    x = []; y = [];
    for i=1:360
        % if i >= sp_t(2) && i <= sp_t(3) % should be || if sp_angle goes through 0
        if sp_range(i) == 1
            m_angle = find(sp_angle == i);
            if size(m_angle,2) == 0,
                continue,
            end
            x = [x; min(sp_xsp(m_angle))]; % was: mean
            y = [y; min(sp_ysp(m_angle))];
       else
            m_angle = find(b_angle == i);
            [m_min, m_pos] = min(border(m_angle, 6));
            m_pos = m_pos + min(m_angle) - 1;
            x = [x; border(m_pos,1)];
            y = [y; border(m_pos,2)];
        end        
    end

    if verbose > 1
        subplot(3,1,2);
        hold on;
        colormap(gray);
        image(img2t); 
        % plot(x,y,'o-');
        plot(x,y,'-');
        hold off;
    end
    
else
    % no touching embryos - simply take outline
    if verbose > 0, fprintf('Detected single embryo, calculating outline ...\n'); end;
    score = 1;
    
    x = []; y = [];
    for i=1:360
        m_angle = find(b_angle == i);
        [m_min, m_pos] = min(border(m_angle, 6));
        m_pos = m_pos + min(m_angle) - 1;
        x = [x; border(m_pos,1)];
        y = [y; border(m_pos,2)];
    end

end

if verbose > 0, fprintf('Retesting for multiple embryos ... '); end;

% find if there are multiple embryos
bord_px = zeros(4,1);
bord_px(1) = length(find(y < 3));
bord_px(2) = length(find(y > nr - 2));
bord_px(3) = length(find(x < 3));
bord_px(4) = length(find(x > nc - 2));

membryo  = 0;
if max(bord_px) >= Pbord_px
    if verbose > 0, fprintf('found\n'); end;
else
    if verbose > 0, fprintf('not found - should be OK\n'); end;
    score = bitset(score,2); % passed final QC
end


if verbose > 0, fprintf('Refining borders with snakes ...\n'); end;
[x, y] = snake(x, y, Param_SNalpha, Param_SNgamma, img2t, Param_SNiter);
x = fix(x);
y = fix(y);

if verbose > 1
    if s_touch == 1, subplot(3,1,3); end;
    hold on;
    colormap(gray);
    image(img2t); plot(x,y,'-');
    hold off;
end

if verbose > 0, fprintf('Final score = %d\n', score); end;
return;

%----------------------------------------------------------------------
% Snake functions

function [xs, ys] = snake(xs, ys, alpha, gamma, image, N);
% ADJSHOW   Repeatedly adjust a snake.
% [XS, YS] = ADJSHOW(XS, YS, ALPHA, GAMMA, IMAGE, N)
%   repeatedly calls ADJUSTSNAKE and displays the snake after each
%   iteration. Arguments and results are as for ADJUSTSNAKE, with the
%   addition of IMAGE on top of which the snake is drawn, and N, the
%   number of iterations to perform.

[xg yg] = gradients(image,3);

for i=1:N;
    [xs, ys] = adjustsnake(xs, ys, xg, yg, alpha, gamma);
end;

%----------------------------------------------------------------------

function [newxs, newys] = adjustsnake(xs, ys, xf, yf, alpha, gamma);
%ADJUSTSNAKE  Make one set of adjustments to snake's control points
%   [NEWXS, NEWYS] = ADJUSTSNAKE(XS, YS, XF, YF, ALPHA, GAMMA)
%   takes vectors of x- and y-coords, XS and YS, and returns updated
%   vectors having made adjustments based on internal elastic forces and
%   image forces.
%
%   The elastic force moves each point towards the mid-point of its
%   neighbours; the distance moved is ALPHA times the distance to
%   the mid-point.
%
%   The external force moves a point at (XI, YI) a distance along the
%   x-axis of GAMMA*XF(XI,YI) and along the y-axis of GAMMA*YF(XI,YI).
%   XF and YF are 2-D structures holding (typically) the x- and
%   y-gradients of the image.

newxs = xs + adj(xs, xs, ys, xf, alpha, gamma);   % do x-coords
newys = ys + adj(ys, xs, ys, yf, alpha, gamma);   % and y-coords



function j = adj(as, xs, ys, af, alpha, gamma);
% Function to do either the x or the y adjustments.

[N, u] = size(as);                      % get no. of control points

% The shifted vectors below are padded by circulating the coords from
% one end of the snake onto the other.
am1 = [as; as(1:2)];                    % unshifted vector
a0 = [as(N); as; as(1)];                % shifted one place down
ap1 = [as(N-1:N); as];                  % shifted two places down

% The elastic adjustment is calculated using the shifted vectors to
% compute a function of each coord and its nearest neighbours.
adja = alpha * (0.5*(am1 + ap1) - a0);
adja = adja(2:N+1);         % remove the padding

% Next line converts the snake coords to linear indices into the
% image force array.
indices = sub2ind(size(af), round(ys), round(xs));

% The image force adjustment is applied simply
j = adja + gamma * af(indices);

%----------------------------------------------------------------------

function [xg, yg] = gradients(image, sigma);
%GRADIENTS Estimate the spatial grey-level gradients
%   [XG, YG] = GRADIENTS(IMAGE, SIGMA) smooths the image with Gaussian
%   of width SIGMA, then convolves it with horizontal and vertical
%   centred differencing masks to get XG and YG respectively.

% Gaussian mask. 6*sigma is quite a generous size which ensures the
% values at the edge of the mask are small
smask = fspecial('gaussian', ceil(6*sigma), sigma);

% Use filter2 rather than conv2 for this convolution because the mask
% is separable and filter2 takes advantage of this to compute faster.
% Use 'same' to keep the output same size as input.
smth = filter2(smask, image, 'same');

% Hor and vert masks - use [1 0 -1]/2 to get the average gradient
% over 2 pixels
xg = conv2(smth, [0.5, 0, -0.5], 'same');
yg = conv2(smth, [0.5; 0; -0.5], 'same');
