% Extract the pixels in the mask (here, "1" denotes ROI pixels, while
% "0" denotes background pixels).
function imV=extract(im,mask)
    im=reshape(im,[],size(im,3));
    imV=im(find(mask),:);
end
