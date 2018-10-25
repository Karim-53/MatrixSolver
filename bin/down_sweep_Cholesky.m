function [y]=down_sweep_Cholesky(A,x)
%down sweep :R�soudre le syst�me triangulaire inf�rieur
[m,n]=size(A);
if (m~=n)
disp('Matrice n''est pas carr�e');
else
y=x;
y(1)=y(1)/A(1,1);
for (i=2:n)
y(i)=y(i)-A(i,1:i-1)*y(1:i-1);
y(i)=y(i)/A(i,i)
end
end