function X=MOL2MASS_CONV3(Y)

MM=MolarMassGERG3(Y);
Y=Y/100;
% [a b c]=size(x);
load('MMiGERG.mat','MMiGERG');
% global MMiGERG;
MMiGERG=MMiGERG';
%Y=permute(Y,[2,1,3]);

dimn=size(Y,1);
for ii=1:size(Y,3)
       
    X(:,:,ii)=reshape(Y(:,:,ii),dimn,21).*MMiGERG./MM(ii,:)'*100;
    
   end