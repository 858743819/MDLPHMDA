function KM=similarity_microbe(interaction)
 
  [~,nm]=size(interaction);
   for i=1:nm
sm(i)=norm(interaction(:,i))^2;
   end
    gamam=nm/sum(sm');
    for i=1:nm
        for j=1:nm
   KM(i,j)=exp(-gamam*(norm(interaction(:,i)-interaction(:,j)))^2);        %calculate Gaussian kernel for the similarity between miRNA: km
       end
    end 
    
