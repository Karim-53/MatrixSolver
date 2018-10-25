function qr_val_propre();
clear all
close all
clc

disp('A initial:')
% A = [7 8 9;4 5 6;1 2 3]
A = magic(4)
disp('valeurs propres avec la fonction matlab:')
matval = sort(eig(A),'ascend')


epsilon = 1e-5;

uvp =  uniquetol (abs(matval) , epsilon  )

if length( uvp  )< length(matval)
    warning('On ne peut pas continuer car 2 valeurs propres ont le meme module')
end


myval = VP(A,epsilon);

disp ([ 'epsilon =' num2str(epsilon) ])
disp('_')
disp('valeurs propres avec la methode QR:')
myval

end

function [vp] = VP(An,epsilon)
epsilon = max(eps,epsilon);
Nmax = 53%nombre max d'itterations


mx = epsilon + 1;
nb=0;
while (mx>epsilon)
    
%      subplot(7,2,nb+1)
%     bar3h(log1p(An(:,end:-1:1)))
%     title([ 'A' num2str(nb) ])
    
    [Q,R]=qr(An);
    An=R*Q;
    mx = elmt_max(An);
    nb=nb+1;
    
    
if (nb>Nmax-1)
    warning(['Arret apres ' num2str(Nmax) ' itterations ##################' ])
    mx=0;
end
end
disp(['nombre d itteration : ' num2str(nb) ])
disp('dernier An calculé :')
An
vp = sort(diag(An),'ascend');
end

function mx=elmt_max(A)
v = [];%initiation
for (i=1:size(A,1)-1)
    v =   [ v ; A(i+1:end,i) ];
end
mx = max( abs(v) );
end