clc
load fisheriris
X       = meas;
[N,D]   = size(X);
[~,~,y] = unique(species);
Xc      = X-repmat(mean(X),N,1);
Xn      = Xc./repmat(max(abs(X)),N,1);
d       = 2;

% PCA CONVENCIONAL
A_euclidean = Xc*Xc';
[V,D]       = eig(1/(N-1)*A_euclidean);     % espectral analysis    
[~,ind]     = sort(diag(D),'descend');      
Y           = V(:,ind(1:d));                % choose maximun eigenvectors
Y           = Y./repmat(max(abs(Y)),N,1);   % normalize between 0 and 1

% PCA MAHALANOBIS
Sigma=cov(X);
A_mahalano =X*(inv(Sigma))*X';%prueba con inversa y sin inversa
[V2,D2]    =eig(A_mahalano);
[~,ind2]   = sort(diag(D2),'descend');
Y2         = V2(:,ind2(1:d));    %In this case the spectral analysis of X*inv(Sigma)*X' already  contains the proyected data
Y2         = Y2-repmat(mean(Y2),N,1);
Y2         = Y2./repmat(max(abs(Y2)),N,1);
%Ploting results from both PCAs

subplot(121)
scatter(Y(:,1),Y(:,2),50,y)
subplot(122)
scatter(Y2(:,1),Y2(:,2),50,y)
