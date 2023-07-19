function Y=MASS2MOL_CONV3(X,MM)

%MM=MolarMassGERG1(Y);
X=X/100;
% [a b c]=size(x);
load('MMiGERG.mat','MMiGERG');
MMiGERG=MMiGERG';
%Y=permute(Y,[2,1,3]);

%X=Y.*MMiGERG./MM'*100;

dimn=size(X,1);
for ii=1:size(X,3)
       
%     X(:,:,ii)=reshape(Y(:,:,ii),dimn,21).*MMiGERG./MM(ii,:)'*100;
    Y(:,:,ii)=reshape(X(:,:,ii),dimn,21)./MMiGERG.*MM(ii,:)'*100;
end
end

% Y=X./MMiGERG.*MM'*100;
