function [A]= Cholesky_fact(A)
%Factorisation Cholesky
[m,n]=size(A);
for (i=1:n)
if (A(i,i)<=0)
disp('Pivot négatif')
else
A(i,i)=sqrt(A(i,i));
A(i+1:m,i)=A(i+1:m,i)/A(i,i);
for (j=i+1:n)
A(j:m,j)=A(j:m,j)-A(j:m,i)*A(j,i);
end ;
end ;
end ;