function [X,E] = lrr(Adaj_dm,alpha1)
J=zeros(292,292);
X=zeros(292,292);
E=zeros(39,292);
Y1=zeros(39,292);
Y2=zeros(292,292);
max_mu=10e10;
mu=10e-4;
rho=1.1;
epsilon=10e-6;
stopC=1;
ATA=Adaj_dm'*Adaj_dm;
inv_x=inv(ATA+eye(292));
while stopC>epsilon
  [U,sigma,V]=svd(X+Y2/mu,'econ');
  sigma=diag(sigma);
  svp=length(find(sigma>1/mu));
  if svp>=1
      sigma=sigma(1:svp)-1/mu;
  else
      svp=1;
      sigma=0;
  end
      J=U(:,1:svp)*diag(sigma)*V(:,1:svp)';
      X=inv_x*(ATA-Adaj_dm'*E+J+(Adaj_dm'*Y1-Y2)/mu);
       temp1 = Adaj_dm-Adaj_dm*X;
       E = solve_l1l2(temp1+Y1/mu,alpha1/mu);
      stopC=max(max(max(abs(temp1-E))),max(max(abs(X-J))));
       Y1=Y1+mu*(temp1-E);
      Y2=Y2+mu*(X-J);
      mu=min(max_mu,rho*mu);
  
end

