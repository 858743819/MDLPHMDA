function KD2=similarity_disease(interaction)

 [nd,~]=size(interaction);
for i=1:nd
sd(i)=norm(interaction(i,:))^2;                                              %calculate gamad for Gaussian kernel calculation
   end
    gamad=nd/sum(sd');
    for i=1:nd
        for j=1:nd
    KD2(i,j)=exp(-gamad*(norm(interaction(i,:)-interaction(j,:)))^2);        %calculate gamad for Gaussian kernel calculation
       end
    end
                

%      for i=1:nd
%        for j=1:nd
%           KD(i,j)=V1(i,j)./sqrt(sum(V1(i,:))*sum(V1(j,:)));
% end
%      end
%       
   
