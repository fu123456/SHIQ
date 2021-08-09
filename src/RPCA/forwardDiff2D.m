function [DhF DvF] = forwardDiff2D(F)
% Compute Forward Finite Difference in 2 D, say DF (Df)
    DhF = diff(F,1,2);
    DhF(:,end+1,:) = 0;
    DvF = diff(F,1,1);
    DvF(end+1,:,:) = 0;
end