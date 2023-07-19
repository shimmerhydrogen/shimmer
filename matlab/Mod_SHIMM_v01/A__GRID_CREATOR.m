clear all
clc
close all
% 
INPUT_b_inner=xlsread('INPUT_Pambour1','PIPE');
INPUT_n_inner=xlsread('INPUT_Pambour1','PIPE_N','B3:Z5');

% TIME=0.25*3600;%s
TIME=2*3600;
dt=1*60;
tt=[0:dt:TIME];
dimt=length(tt);


Ainput_inner=INPUT_b_inner(:,1:3);% |--n°Branch--|--NodeIN--|--NodeOUT--|
LL_inner=INPUT_b_inner(:,4)*1e3;%m
DD_inner=INPUT_b_inner(:,5); %m
epsi_inner=INPUT_b_inner(:,6)*1e-3;%m
COMP_inner=INPUT_b_inner(:,7);
REG_inner=INPUT_b_inner(:,8);
VAL_inner=INPUT_b_inner(:,9);
RES_inner=INPUT_b_inner(:,10);
nGridPoints=INPUT_b_inner(:,11);

HH_inner=INPUT_n_inner(:,2);%m

%  dimb0=size(INPUT_b_inner,1);
%  dimn0=max(max(INPUT_b_inner(:,2:3)));
% 
% BR0 = INPUT_b_inner(:,1);
% IN0 =  INPUT_b_inner(:,2); 
% OUT0 = INPUT_b_inner(:,3);
% Asp_p0=sparse(IN0,BR0,ones(1,dimb0),dimn0,dimb0);
% Asp_m0=sparse(OUT0,BR0,-ones(1,dimb0),dimn0,dimb0);
% Asp0=Asp_p0+Asp_m0;
% 
% delH_inner=-Asp0'*HH_inner;

IN=[];
OUT=[];
CORR=[];
nodi_notevoli=[];
nodi_notevoli_corr=[];
INPUT_b=[];
% INPUT_n=[];
COMP=[];
REG=[];
VAL=[];
RES=[];
NP=[];
PIPE=[];
in=0;
for kk=1:length(INPUT_b_inner(:,1))
    
    branch(kk).xx(1,:)=[0:LL_inner(kk)./nGridPoints(kk):LL_inner(kk)];
    %nodi notevoli
    branch(kk).xx(2,1)=INPUT_b_inner(kk,2); %in
    branch(kk).xx(2,end)=INPUT_b_inner(kk,3); %out
    %nodi
    branch(kk).xx(3,:)=[in+1:in+length(branch(kk).xx)];
        
    if kk>=2
       if any(branch(kk).xx(2,1)==nodi_notevoli) %ingresso
           pos=find(branch(kk).xx(2,1)==nodi_notevoli);
           branch(kk).xx(3,1)=nodi_notevoli_corr(pos(1));
           pos=[];
           flagIN=1;
       end
       
        if any(branch(kk).xx(2,end)==nodi_notevoli) %uscita
           pos=find(branch(kk).xx(2,end)==nodi_notevoli);
           branch(kk).xx(3,end)=nodi_notevoli_corr(pos(end));
           pos=[];
           flagOUT=1;
        end
        
        if flagIN==1 & flagOUT==1
           branch(kk).xx(3,2:end-1)=[in+1:in+size(branch(kk).xx,2)-2];
           in=branch(kk).xx(3,end-1);
        elseif flagIN==1 & flagOUT==0
           branch(kk).xx(3,2:end)=[in+1:in+size(branch(kk).xx,2)-1];
           in=branch(kk).xx(3,end);
        elseif flagIN==0 & flagOUT==1
           branch(kk).xx(3,1:end-1)=[in+1:in+size(branch(kk).xx,2)-1]; 
           in=branch(kk).xx(3,end-1);
            else
           branch(kk).xx(3,1:end-1)=[in+1:in+size(branch(kk).xx,2)-1];
           in=branch(kk).xx(3,end);
        end   
    else
        in=branch(kk).xx(3,end);  
    
    end
            
    
    
    nodi_notevoli=[nodi_notevoli;branch(kk).xx(2,:)'];
    nodi_notevoli_corr=[nodi_notevoli_corr;branch(kk).xx(3,:)'];
    NN=[nodi_notevoli,nodi_notevoli_corr];
    
    flagIN=0;
    flagOUT=0;
            
%     NODE_BR(kk).IN=branch(kk).xx(3,1:end-1);
%     NODE_BR(kk).OUT=branch(kk).xx(3,2:end);
    IN=[IN;branch(kk).xx(3,1:end-1)']; 
    OUT=[OUT;branch(kk).xx(3,2:end)'];

    
    LL_pipe(kk)=LL_inner(kk)./nGridPoints(kk);
    lung(kk)=length(branch(kk).xx);
    branch(kk).xx(4,:)=LL_pipe(kk);
    branch(kk).xx(5,:)=DD_inner(kk);
    branch(kk).xx(6,:)=epsi_inner(kk);

    branch(kk).xx(7,:)=linspace(HH_inner(INPUT_b_inner(kk,2)),HH_inner(INPUT_b_inner(kk,3)),nGridPoints(kk)+1);

    altitude(branch(kk).xx(3,:))=branch(kk).xx(7,:);
    
   INPUT_b=[INPUT_b;branch(kk).xx(3,1:end-1)',branch(kk).xx(3,2:end)',branch(kk).xx(4:6,2:end)'];
   
   if COMP_inner(kk)==1
   COMP=[COMP;branch(kk).xx(3,1:end-1)'];
   NP=[NP;branch(kk).xx(3,1:end-1)'];
   elseif REG_inner(kk)==1
   REG=[REG;branch(kk).xx(3,1:end-1)'];
   NP=[NP;branch(kk).xx(3,1:end-1)'];
   elseif VAL_inner(kk)==1
   VAL=[VAL;branch(kk).xx(3,1:end-1)'];
   NP=[NP;branch(kk).xx(3,1:end-1)'];
   elseif RES_inner(kk)==1
   RES=[RES;branch(kk).xx(3,1:end-1)'];
   NP=[NP;branch(kk).xx(3,1:end-1)'];
   else
   PIPE=[PIPE;branch(kk).xx(3,1:end-1)'];
   end
   
   xx_old(kk).xx=branch(kk).xx(1,1:end-1);
    
end
 dimb=size(INPUT_b,1);
 BR=[1:1:dimb]';
 INPUT_b=[BR,INPUT_b];
 dimn=max(max(INPUT_b(:,2:3)));
 
INPUT_n=zeros(dimn,25);
corrispondenze=unique(NN(find(NN(:,1)~=0),:),'rows');
INPUT_n(:,1)=[1:1:dimn];
INPUT_n(corrispondenze(:,2),3:end)=INPUT_n_inner(corrispondenze(:,1),3:end);
INPUT_n(:,2)=altitude;

Ainput=INPUT_b(:,1:3);

BR = INPUT_b(:,1);
IN =  INPUT_b(:,2); 
OUT = INPUT_b(:,3);
Asp_p=sparse(IN,BR,ones(1,dimb),dimn,dimb);
Asp_m=sparse(OUT,BR,[-ones(1,dimb)],dimn,dimb);
Asp=Asp_p+Asp_m;

figure(100)
GRAPH=digraph(IN,OUT);
plot(GRAPH,'layout','force');

GRAPH1=graph(IN,OUT);
gradi=degree(GRAPH1);
INNER_H=find(gradi>2);


save('DATA_INPUT_Trial2023_3.mat')



