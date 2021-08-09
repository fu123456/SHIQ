% Xiaojie Guo, Oct 2013. 
% Questions? xj.max.guo@gmail.com
%
% Reference: Robust Separation of Reflection Using Multiple Images,
% Xiaojie Guo, Xiaochun Cao and Yi Ma, CVPR 2014
function [Fot, T, R, N,  transformations] =  SID(I0, transformations, canonSize, mode)
% ---------------------------------------------------
% Superimposed Image Decomposition
% \|L\|_* + lambda1*\|M\|_1 + lambda2*\|N\|_1+\lambda3*\|K\|_1 +
% \lambda4*\|Q\|_1 + \lambda5\|K@Q\|_1 + \lambda6\|Omega-K-Q\|_F^2
% s.t. Fot = T+M; M = N+R; L=T; K = DT; Q = DR; T>=0; R>=0;
% ---------------------------------------------------
if nargin < 4
    mode = 1;
end
% parametric tranformation model
para.transformType = 'HOMOGRAPHY'; 
% one of 'TRANSLATION', 'EUCLIDEAN', 'SIMILARITY', 'AFFINE','HOMOGRAPHY'
% main loop
para.stoppingDelta = .1; % stopping condition of main loop
if mode == 1
    para.maxIter = 10; % maximum iteration number of main loops
else
    para.maxIter = 1; 
end

% inner loop
nbOfFrames = length(transformations);
%----------------------------------------
% Here, lambda1 to lambda6 are the coefficients can be tuned
%----------------------------------------
dv = canonSize(1);
dh = canonSize(2);
coef.lambda1 = 0.3/sqrt(dv*dh);
coef.lambda2 = 50/sqrt(dv*dh);
coef.lambda3 = 1/sqrt(dv*dh);
coef.lambda4 = 5/sqrt(dv*dh);
coef.lambda5 = 50/sqrt(dv*dh);
coef.lambda6 = 50/sqrt(dv*dh);
coef.imgSize = canonSize;
%----------------------------------------
I0x = cell(1,nbOfFrames);
I0y = I0x;
for fileIndex = 1 : nbOfFrames
    I0x{1,fileIndex} = imfilter( I0{1,fileIndex}, (-fspecial('sobel')') / 8 );
    I0y{1,fileIndex} = imfilter( I0{1,fileIndex},  -fspecial('sobel')   / 8 ); 
end

%% get the initial input images in canonical frame
xi_initial = cell(1,nbOfFrames) ; % initial transformation parameters
for i = 1 : nbOfFrames
    xi_initial{i} = projective_matrix_to_parameters('HOMOGRAPHY',transformations{i});
end

%% start the main loop
iterNum = 0 ;  % iteration number of outer loop
converged = 0 ;
prevObj = inf ; % previous objective function value 
xi = cell(1,nbOfFrames) ;
while ~converged
    iterNum = iterNum + 1 ;
    Fot = zeros(dv*dh, nbOfFrames);
    J = cell(1,nbOfFrames);
    disp(['Outer Loop Iter ' num2str(iterNum)]) ;
    for fileIndex = 1 : nbOfFrames
        Tfm = fliptform(maketform('projective',transformations{fileIndex}'));
        I   = vec(imtransform(I0{1,fileIndex}, Tfm,'bicubic','XData',...
            [1 canonSize(2)],'YData',[1 canonSize(1)],'Size',canonSize));
        Iu  = vec(imtransform(I0x{1,fileIndex},Tfm,'bicubic','XData',...
            [1 canonSize(2)],'YData',[1 canonSize(1)],'Size',canonSize));
        Iv  = vec(imtransform(I0y{1,fileIndex},Tfm,'bicubic','XData',...
            [1 canonSize(2)],'YData',[1 canonSize(1)],'Size',canonSize));
        y   = I; %vec(I);
        Fot(:,fileIndex) = y ;
        % transformation matrix to parameters
        xi{fileIndex} = projective_matrix_to_parameters...
            (para.transformType,transformations{fileIndex}) ; 
        % Compute Jacobian
        J{1,fileIndex} = image_Jaco(Iu, Iv, canonSize, ...
            para.transformType, xi{fileIndex});
    end

       
    % Superimposed Image Decomposition inner loop
    % -----------------------------------------------------------------
    % -----------------------------------------------------------------
    % using QR to orthogonalize the Jacobian matrix
    QR_Q = cell(1,nbOfFrames);
    QR_R = QR_Q;
    for fileIndex = 1 : nbOfFrames
        [QR_Q{fileIndex}, QR_R{fileIndex}] = qr(J{1,fileIndex},0);
    end
    [L,M,N,K,Q,T,R,delta_xi,Omega] = SID_inner(Fot,QR_Q,coef);
    
    for fileIndex = 1 : nbOfFrames
        delta_xi{fileIndex} = (QR_R{fileIndex})\delta_xi{fileIndex} ;
    end
   % step in paramters
    for i = 1 : nbOfFrames
        xi{i} = xi{i} + delta_xi{i};
        transformations{i} = parameters_to_projective_matrix(para.transformType,xi{i});
    end
    % -----------------------------------------------------------------
    curObj = norm(svd(L),1) + coef.lambda1*sum(abs(M(:))) + coef.lambda2*sum(abs(N(:)))...
        + coef.lambda3*sum(abs(K{1}(:))+abs(K{2}(:))) + coef.lambda4*sum(abs(Q{1}(:))+abs(Q{2}(:)))...
        + coef.lambda5*sum(abs(K{1}(:).*Q{1}(:))+abs(K{2}(:).*Q{2}(:))) +...
        + coef.lambda6*norm(Omega{1} - K{1} - Q{1},'fro')^2+coef.lambda6*norm(Omega{2} - K{2} - Q{2},'fro')^2;
    disp(['Previous objective function: ' num2str(prevObj) ]);
    disp(['Current objective function : ' num2str(curObj) ]);


    if curObj>prevObj
        L = L_old;
        M = M_old;
        N = N_old;
        K = K_old;
        Q = Q_old;
        T = T_old;
        R = R_old;
        transformations = transformations_old;
    else
        L_old = L;
        M_old = M;
        N_old = N;
        K_old = K;
        Q_old = Q;
        T_old = T;
        R_old = R;
        transformations_old = transformations;
    end
    if ( (prevObj - curObj < para.stoppingDelta) || iterNum >= para.maxIter )
        converged = 1;
        if ( prevObj - curObj >= para.stoppingDelta )
            disp('Maximum iterations reached') ;
        end
    else
        prevObj = curObj;
    end

end


