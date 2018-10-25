function [x]=remontee(A,b)
n=size(A,1);
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
    x(i)=b(i)/A(i,i);
    for j=n:-1:i+1
        x(i)=x(i)-x(j)*A(i,j)/A(i,i);
   end
end
x=transpose(x);

end

