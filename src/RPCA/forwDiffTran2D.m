function [DtQ]= forwDiffTran2D(Qh, Qv)
% Compute Transpose of FFD, say D'Q
    [h w d] = size(Qh);
    DtQ = cat(2, zeros(h,1,d), -diff(Qh,1,2))+...
        cat(1, zeros(1,w,d), -diff(Qv,1,1));%DhtQ + DvtQ;
end