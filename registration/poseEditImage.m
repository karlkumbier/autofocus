function out_image = poseEditImage(out_image,inverted)
% inverted: 1 ap_inverted; 2: dv_inverted; 3: ap_dv_inverted
if inverted == 1
        out_image(:,:,1) = fliplr(out_image(:,:,1));
        out_image(:,:,2) = fliplr(out_image(:,:,2));
        out_image(:,:,3) = fliplr(out_image(:,:,3));
elseif inverted == 2
        out_image(:,:,1) = flipud(out_image(:,:,1));
        out_image(:,:,2) = flipud(out_image(:,:,2));
        out_image(:,:,3) = flipud(out_image(:,:,3));
elseif inverted == 3
        out_image(:,:,1) = fliplr(out_image(:,:,1));
        out_image(:,:,2) = fliplr(out_image(:,:,2));
        out_image(:,:,3) = fliplr(out_image(:,:,3));
        out_image(:,:,1) = flipud(out_image(:,:,1));
        out_image(:,:,2) = flipud(out_image(:,:,2));
        out_image(:,:,3) = flipud(out_image(:,:,3));
end
