% Using the multi-time Otsu algorithm to produce highlight mask, since
% the only-one-time one often generates highlight masks with much
% noise even errors
function [mask] = Rolling_L_Otsu(image, iter_num)
    Lab=rgb2lab(image);
    L=Lab(:,:,1);
    L=(L-min(L(:)))./(max(L(:))-min(L(:)));
    mask=ones(size(image,1),size(image,2));
    for i=1:iter_num
        L_temp=extract(L,mask);
        T=graythresh(L_temp);
        mask(find(L<=T))=0;
    end
end
