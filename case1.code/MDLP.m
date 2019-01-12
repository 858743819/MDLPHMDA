 A=importdata('disease-microbe associationg number.xlsx');                                     %A: Binary relations between disease and microbe, 1st column:microbe, 2nd column:disease
 A=A.Sheet1;
 mi_num=importdata('microbe identifier.xlsx');
 d_num=importdata('disease identifier.xlsx');
  KKD=importdata('disease symptom similarity.xlsx','Sheet1');  
  alpha=0.3;                                                                % define the parameters
  maxit=10e-06;                                                             
  nmic=max(A(:,2));                                                         % nd:the number of diseases
  nd=max(A(:,1));                                                           % nm:the number of microbes
  [pp,qq]=size(A);                                                          % pp:the number of known diseae-microbe associations
  adaj_dm=zeros(nd,nmic);   
  alpha1=0.1;
  KD1=zeros(nd,nd);
   predict_di_mi_score=cell(1,nd);
  uu=cell(1,3);
  allranking=[];
for i=1:61
    KD1(KKD(i,1),KKD(i,2))=KKD(i,3);  
    KD1(KKD(i,2),KKD(i,1))=KKD(i,3);                                           %adaj_dm: adajency matrix for the disease-microbe association network
end
                                 
for i=1:pp
        adaj_dm(A(i,1),A(i,2))=1;                                           %adaj_dm: adajency matrix for the disease-miRNA association network
end

    
            Adaj_dm=adaj_dm;                                                %assign adaj_dm to Adaj_dm
            Adaj_dm(A(i,1),A(i,2))=0;                                      %a known association is regarded as an unknown, the corresponding matrix of the 1 is changed to 0
            AAdaj_dm=Adaj_dm;
            [X,E] = lrr(AAdaj_dm,alpha1);                                   %SLM
            AAdaj_dm= AAdaj_dm*X;   
            ADAJ_DM= AAdaj_dm;
            ADAJ_sco= AAdaj_dm;
            ADAJ_M= AAdaj_dm';
            KM=similarity_microbe(  AAdaj_dm); 
            KD2=similarity_disease(AAdaj_dm);
            KD= (KD1+KD2)./2;     
            change=1;       
            Norm_adaj_dm1=KD;  
            Norm_adaj_m1=KM;
            d1=sum( Norm_adaj_dm1,2);                                 
            d2=sum( Norm_adaj_m1,2);
            D1=diag(1./sqrt(d1));
            D1(D1==Inf)=0;                                                  
            D2=diag(1./sqrt(d2));
            D2(D2==Inf)=0;
            Norm_adaj_dm2=D1* Norm_adaj_dm1*D1;                              % normalization the integrated similarity of KD
            Norm_adaj_m2=D2* Norm_adaj_m1*D2;                                % normalization the integrated similarity of KM  
        while( change>maxit)    
           
       NEW_adaj_dm=alpha* Norm_adaj_dm2* ADAJ_DM+(1-alpha)* Adaj_dm;          %the matrix of predicting score between disease-miRNA
       NEW_adaj_m=alpha* Norm_adaj_m2* ADAJ_M+(1-alpha)* Adaj_dm';            %the matrix of predicting score between miRNA-lncRNA  
       New_adaj_sco=0.5*NEW_adaj_dm+0.5*NEW_adaj_m';                         %predicted score of miRNA-disease associations
       change=sum(sum((abs(New_adaj_sco- ADAJ_sco))));                       % the iterative change     
       ADAJ_DM=NEW_adaj_dm;                                                  %assign NEW_adaj_dm to ADAJ_DM
       ADAJ_M= NEW_adaj_m;                                                   %assign NEW_adaj_m to ADAJ_M
       ADAJ_sco=0.5.*ADAJ_DM+0.5.*ADAJ_M';
        end
        %golbal
         Sco= ADAJ_sco;                                                         %the matrix of predicting score between disease-miRNA
    for oo=1 :nd
  disease_rna_num=adaj_dm(oo,:);
  disease_rna_score= Sco(oo,:);
  irrevelant_pair_num=find(disease_rna_num==0);                              
  irrevelant_pair_score=disease_rna_score(irrevelant_pair_num);                 %row vector
  [er,irrevelant_num]=size(irrevelant_pair_score);
   pair=zeros(irrevelant_num,3);
   pair(:,1)=oo;
   pair(:,2)=irrevelant_pair_num';
   pair(:,3)=irrevelant_pair_score';
   di_rna_score=sortrows(pair,-3);
  each_disease=cell(irrevelant_num+1,3);
  each_disease{1,1}='Disease Name';
  each_disease{1,2}='microbe Name';
  each_disease{1,3}='Score';
   each_disease2=cell(irrevelant_num,3);
        for j=1:irrevelant_num
            each_disease{j+1,1}=d_num.textdata.Sheet1(di_rna_score(j,1));           
            each_disease{j+1,2}=mi_num.textdata.Sheet1(di_rna_score(j,2));          
            each_disease{j+1,3}=di_rna_score(j,3);                           
        end
        predict_di_mi_score{1,oo}=each_disease;                             % predict_di_mi_score is the scores between candidate micirbes and disesases
       

       for g=1:irrevelant_num
            each_disease2(g,1)=d_num.textdata.Sheet1(di_rna_score(g,1));           
            each_disease2(g,2)=mi_num.textdata.Sheet1(di_rna_score(g,2));          
            each_disease2{g,3}=di_rna_score(g,3);     
               
        end
   
     predict_di_mi_score2{1,oo}=each_disease2;
     uu=[uu;predict_di_mi_score2{1,oo}];                                    %uu is the socres of all potential microbe-disease pairs
   
    end
 xlswrite('all.xlsx',uu);

for i=1:nd
      d_char_name=d_num.textdata.Sheet1(i);                                  %return ith disease name;class:cell, use char()
      str=['./result\',char(d_char_name),'.xlsx'];                          
      d_char_data=[predict_di_mi_score{1,i}{:}];
      d_char_data=reshape(d_char_data,length(d_char_data)/3,3);
      xlswrite(str,d_char_data);
end

