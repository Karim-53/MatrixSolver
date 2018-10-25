%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parametres  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[]=main_4;
clc;
clf
close('all');
%fonction
f = inline('2*sin(x+y)','x','y')
%dimension
n = 9

% intervalle ]0 l[
l = 2*pi


[b h X Y]= vec_b(f,n,l);

A = mat_A( n, h )


[L, U] = LUdecompCrout(A)

%Solve Lv=c par une descente
v=descente (L,b)


%Solve Uu=v par une remontée
u=remontee(U,v)




uu = A\b;
% uu
% u
%max(uu-u,1e-9)-1e-9
disp('erreur relative (comparé au programme matlab): ')
norm(u-uu)
%TODO
%erreur relative theorique



%affichage
figure('Name','Résultat');
mesh(X,Y(end:-1:1,1:end),To_Mat(u))

disp('end.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Drive Test  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pour tester la fonction qui genere la matrice A (main drive)
function[]=Drive_mat_A;
for(n=2:6)
A =  mat_A( n );
if (  (n-1)^2 ~= length(A)  )
   disp('erreur: la taille de la matrice A est diff de (n-1)²');
end

end
disp('(end test A)')
end

function[]=Drive_decomp;
A = [[2 4 4]',[1 3 1]',[1 5 6]']'
LL = [[1.0 0 0]',[0.5 1.0 0]',[0.5 3.0 1.0]']';
UU = [[2 4 4]',[0 1 -1]',[0 0 7]']';

[L U]=LUdecompCrout(A)

if LL~=L
    disp('erreur l');
    L
    LL
end
if UU~=U
    disp('erreur u');
    U
    UU
end
if ( A~=L*U )
    disp('A est diff de  L*U' )
end
end
function[]=Drive_To_vector;
b = [ 7 8 9;4 5 6;1 2 3 ];
u=  To_vec(b);
if u'~=1:9
    u
    disp('erreur 1')
    pause
else
    disp('OK')
end
bb =  To_Mat(u);
if bb~=b
    bb
    disp('erreur 2')
    pause
else
    disp('OK')
end







end


function[]=Drive_resolution;
A = [3 6;5 7]'
b = [9,4]'
xx = A\b

[L U]=LUdecompCrout(A);

%Solve Lv=c par une dessente
v = L\(b);

%Solve Uu=v par une remontée
u = U\v
if xx~=u
    xx-u;
    disp('erreur')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ A ] = mat_A( n , h)
if nargin < 2
  h = 1%default
end
B=diag( ones(n-2,1),-1 )*(-1) + diag( ones(n-2,1),1 )*(-1) + diag( ones(n-1,1) )*4;
tmp = ones(0,0);%matrice vide
for( i=1:n-1)%Construction de A par bloc
tmp = blkdiag(tmp,B);
end
%Ajout des identité
A = tmp+ diag( ones( (n-1)*(n-2) ,1),   n-1 )*(-1) ;
A = A  + diag( ones( (n-1)*(n-2) ,1), -(n-1) )*(-1) ;
A = A./(h*h);




if (det(A)==0)
   disp('la matrice A est non inversible. on peut pas appliquer la methode de LU')
    return;
end
if ( size(A,1)~=size(A,2) )
   disp('la matrice A n`est pas carré')
   return;
end






%TODO
%A=spalloc(N,N,5*N); % matrice creuse de 5 elts maxi / ligne
end


function [b h X Y]=  vec_b(f,n,L)
h=L/n
i= 1:(n-1);
X0 = h*i;


b = f(X0, flip(X0)' );


figure('Name','b sous forme matriciel (n=33)');
[X,Y] = meshgrid(X0);
mesh(X,Y(end:-1:1,1:end),b)
title([ 'f(x,y)= '  formula(f) ]);

%valeurs au limites
% un = f(X0,L )%nord
% us = f(X0,0 )%sud
% 
% ue = f(L,flip(X0)' )%est
% uo = f(0,flip(X0)' )%ouest
% 
% 
% b(1  ,1:n-1) =b(1  ,1:n-1) + un./h^2;
% b(n-1,1:n-1) =b(n-1,1:n-1) + us./h^2;
% 
% b(1:n-1,1  ) =b(1:n-1,1  ) + uo./h^2;
% b(1:n-1,n-1) =b(1:n-1,n-1) + ue./h^2;

b =  To_vec(b);
end


function [b]=  To_vec(b)
    b = flip( reshape( b(1:end,end:-1:1)' ,[],1) );
end
function [d]=  To_Mat(u)
    c=reshape( flip(u) ,sqrt(size(u,1)),sqrt(size(u,1)) )';
    d = c(1:end, end:-1:1);
end



function [L, U] = LUdecompCrout(A)
m=size(A,1); % pour A carré
L=eye(m);%diag à 1
U=zeros(m);

U(1,:)=A(1,:);
L(:,1) = A(:,1) / U(1,1);

for j=2:m
   for i=j:m
       R=0;
       for k=1:j 
         R=R+(L(j,k)*U(k,i));
       end
         U(j,i)=A(j,i)-R;
   end 
   for i=j+1:m 
        R=0 ;
        for k=1:j 
            R=R+L(i,k)*U(k,j);
        end 
        L(i,j)=(A(i,j)-R)/U(j,j);
   end 
end
end




function [y]=descente (A,b)
n=size(A,1);
y(1)=b(1)/A(1,1);
for i=1:n
    y(i)=b(i)/A(i,i);
    for j=1:i-1
        y(i)=y(i)-A(i,j)*y(j)/A(i,i);
    end
end
y=y';

end

function [x]=remontee(A,b)
n=size(A,1);
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
    x(i)=b(i)/A(i,i);
    for j=n:-1:i+1
        x(i)=x(i)-x(j)*A(i,j)/A(i,i);
   end
end
x=x';

end

