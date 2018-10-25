function [y]=descente (A,b)
n=size(A,1);
y(1)=b(1)/A(1,1);
for i=1:n
    y(i)=b(i)/A(i,i);
    for j=1:i-1
        y(i)=y(i)-A(i,j)*y(j)/A(i,i);
    end
end
y=transpose(y);

end
