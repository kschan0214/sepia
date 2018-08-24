% %%%%%%%%%%%%%%%%%%%%% 
function y = operator4fpc(trans,m,n,x,inds,OMEGA)

if ~trans
   y = A_dct(x,OMEGA);
else
   y = At_dct(x,OMEGA,n); 
end
 

