function lopt = BDU_opt(A,Q,W,b,Ea,Eb,H)
% function lopt = BDU_opt(A,H,Q,R,x0,P0,y)
%
% Compute the optimal parameter gamma for the BDU filter
%
% input:  A,H,Q,W,b -> matrices of the system
%         x0, P0 -> initial conditions      
%         Ea,Eb,H -> uncertainty parameters
%        
% output: lopt -> optimal parameter

%
% Author: Marcos Rogerio Fernandes
% E-mail: eng.marofe@hotmail.com
% Date: 25/10/2018

if (H==0)
    lopt=0;
else
    n=size(H,2);
    l0=1.1*norm(H'*W*H);
    %%
    l=l0:0.01:l0+1;
    for i=1:length(l)
    Ql=Q+l(i)*Ea'*Ea;
    Wl=W+W*H*(l(i)*eye(n)-H'*W*H)\H'*W;
    xl=(Ql+A'*Wl*A)\(A'*Wl*b+l(i)*Ea'*Eb);
    Gl(i)=xl'*Q*xl+l(i)*(Ea*xl-Eb)'*(Ea*xl-Eb)+(A*xl-b)'*Wl*(A*xl-b);
    if i>1
        if Gl(i)>Gl(i-1)
            break
        end
    end
    end
%     plot(l,Gl)
%     hold on
     [Gl0 i0]=min(Gl);
%     plot(l(i0),Gl0,'ro')
lopt=l(i0);
end
end