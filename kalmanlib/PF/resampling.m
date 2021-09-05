function [p,w] = resampling(p0,w)
% from https://www.youtube.com/watch?v=wNQVo6uOgYA
Np=length(w);
p=zeros(size(p0));
i=randi(Np);
Wmax=2*max(w);
beta=0;
for j=1:Np
    beta=beta+Wmax*rand;
    while (w(i)<beta)
        beta=beta-w(i);
        if i<Np
            i=i+1;
        else
            i=1;
        end
    end
    p(j,:)=p0(i,:);
end
w=1/Np*ones(Np,1);
end

