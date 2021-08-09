% Xiaojie Guo, Oct 2013. 
% Questions? xj.max.guo@gmail.com
%
% Reference: Robust Separation of Reflection Using Multiple Images,
% Xiaojie Guo, Xiaochun Cao and Yi Ma, CVPR 2014
function [L,M,N,K,Q,T,R,dt_dual,Omega] = SID_inner(Fot, J, coef)
% Inner Loop of SID
tol = 1e-6;
% More Iterations, More Accurate. Usually, maxIter can be set to 40-100
maxIter = 40;

dt = size(Fot,2);
dv = coef.imgSize(1);
dh = coef.imgSize(2);
sizeF = [dv dh dt];

[Omegah Omegav] = forwardDiff2D(reshape(Fot,sizeF));
Omega = cell(1,2);
Omegah = reshape(Omegah,[dh*dv, dt]);
Omegav = reshape(Omegav,[dh*dv, dt]);
Omega{1} = Omegah;
Omega{2} = Omegav;
T = zeros(dv*dh,dt);
L = T;
M = T;
N = T;
K = cell(1,2);
Kh = T;
Kv = T;
Q = K;
Qh = T;
Qv = T;
R = T;

Z1 = T;
Z2 = T;
Z3 = T;

Z4h = Kh;
Z4v = Kv;

Z5h = Qh;
Z5v = Qv;

dt_dual = cell(1,dt) ;
dt_dual_matrix = zeros(dh*dv, dt) ;

Dh = psf2otf([1,-1],sizeF);
Dv = psf2otf([1;-1],sizeF);
eigDtDs = abs(Dh).^2 + abs(Dv).^2;
F_norm = norm(Fot,'fro');

mu = 0.001;
mu_bar = 1e10;
rho = 1.25;
iter = 0;
converged = false;
imrank = 1;
    s = [];
while ~converged
    iter = iter + 1;
    % Update L
    temp_T = T - Z3*(1/mu);
    [U S V] = svd(temp_T,'econ');
    diagS = diag(S);
    if length(diagS)>imrank;
        diagS(imrank+1:end)=0;
    end
    L = U*diag(pos(diagS-1/mu))*V';
    L(L<0) = 0;
    
    % Update M
    temp_T = 0.5*(Fot+dt_dual_matrix-T+Z1/mu+N+R-Z2/mu);
    M = sign(temp_T).*pos(abs(temp_T)-coef.lambda1/(mu));
    
    % Update N
    N = (Z2+mu*(M-R))/(2*coef.lambda2+mu); 
  
    % Update K
    [Dht Dvt] = forwardDiff2D(reshape(T,sizeF));
    temp_T = (coef.lambda6*(Omegah-Qh)+0.5*mu*reshape(Dht,[dh*dv,dt])-Z4h*0.5)/(coef.lambda6+0.5*mu);
    Kh = sign(temp_T).*pos(abs(temp_T)-(coef.lambda3+coef.lambda5*abs(Qh))/(mu+2*coef.lambda6));
    temp_T = (coef.lambda6*(Omegav-Qv)+0.5*mu*reshape(Dvt,[dh*dv,dt])-Z4v*0.5)/(coef.lambda6+0.5*mu);
    Kv = sign(temp_T).*pos(abs(temp_T)-(coef.lambda3+coef.lambda5*abs(Qv))/(mu+2*coef.lambda6));

    % Update Q
    [Dhr Dvr] = forwardDiff2D(reshape(R,sizeF));
    temp_T = (coef.lambda6*(Omegah-Kh)+0.5*mu*reshape(Dhr,[dh*dv,dt])-Z5h*0.5)/(coef.lambda6+0.5*mu);
    Qh = sign(temp_T).*pos(abs(temp_T)-(coef.lambda4+coef.lambda5*abs(Kh))/(mu+2*coef.lambda6));
    temp_T = (coef.lambda6*(Omegav-Kv)+0.5*mu*reshape(Dvr,[dh*dv,dt])-Z5v*0.5)/(coef.lambda6+0.5*mu);
    Qv = sign(temp_T).*pos(abs(temp_T)-(coef.lambda4+coef.lambda5*abs(Kv))/(mu+2*coef.lambda6));
    
    % Update T
    KZh = Kh+Z4h/mu;
    KZv = Kv+Z4v/mu;
    DtK = forwDiffTran2D(reshape(KZh,sizeF),reshape(KZv,sizeF));   

    
    temp_T = reshape((Fot+dt_dual_matrix-M+L+(Z1+Z3)*(1/mu)),sizeF) + DtK;
    T_tmp = real(ifftn(fftn(temp_T)./(2+eigDtDs)));
    T = reshape(T_tmp,[dv*dh,dt]);
    T(T<0) = 0;
    
    % Update R
    QZh = Qh+Z5h/mu;
    QZv = Qv+Z5v/mu;
    DtQ = forwDiffTran2D(reshape(QZh,sizeF),reshape(QZv,sizeF)); 
    
    temp_T = reshape((M-N+Z2*(1/mu)),sizeF)  + DtQ;
    R_tmp = real(ifftn(fftn(temp_T)./(1+eigDtDs)));
    R = reshape(R_tmp,[dv*dh,dt]);
    R(R<0) = 0;
    
    % Update DeltaGamma
    temp_T = (Fot - T- M +Z1/mu);
    for i = 1 : dt
        dt_dual{i} =  - J{i}'*temp_T(:,i);
        dt_dual_matrix(:, i) = J{i}*dt_dual{i} ;
    end
    
    % Update Largrange Multipliers Z1-Z6
    temp_T = Fot+dt_dual_matrix - T - M;
    Z1 = Z1 + mu*temp_T;
    
    temp_T = M - N- R;
    Z2 = Z2 + mu*temp_T;
    
    temp_T = L- T;
    Z3 = Z3 + mu*temp_T;
    
    [Dht Dvt] = forwardDiff2D(reshape(T,sizeF));
    temp_T = Kh - reshape(Dht,[dh*dv,dt]);
    Z4h = Z4h + mu*temp_T;
    temp_T = Kv - reshape(Dvt,[dh*dv,dt]);
    Z4v = Z4v + mu*temp_T;
    
    [Dhr Dvr] = forwardDiff2D(reshape(R,sizeF));
    temp_T = Qh - reshape(Dhr,[dh*dv,dt]);
    Z5h = Z5h + mu*temp_T;
    temp_T = Qv - reshape(Dvr,[dh*dv,dt]);
    Z5v = Z5v + mu*temp_T;
    
    % Update mu
    mu = min(rho*mu, mu_bar);

    %% stop Criterion    
    stopCriterion = (norm(Fot+dt_dual_matrix-T-N-R,'fro')) / F_norm;
%     s(iter) = stopCriterion;
%     figure(101);clf;plot(1:iter,s,'r-','linewidth',3);
    if mod(iter,10)==0
    disp(['Current Stop Criterion is ', num2str(stopCriterion)])

    figure(1);imshow(reshape((T(:,1)),[dv dh]),[]);
    figure(2);imshow(reshape((R(:,1)),[dv dh]),[]);

    end
    if stopCriterion < tol
        converged = true;
        K{1} = Kh; K{2} = Kv;
        Q{1} = Qh; Q{2} = Qv;
    end    
      
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;      
        K{1} = Kh; K{2} = Kv;
        Q{1} = Qh; Q{2} = Qv;
      
    end
    
end