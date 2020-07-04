function [ys,y,w,b] = saida(x,p,q,s,c,m,n,o) 
%     a = 0;
%     b = 0;
%     y = zeros(m,1);
%     w = ones(m,1);
%     u = zeros(n,m);
%     for j=1:m
%        for i=1+n*(o-1):n*o
%           y(j) = (p(i,j)*x(i-n*(o-1))) + q(j) + y(j);
          st = 1+n*(o-1):n*o;
%           y(j) = sum(p(1+n*(o-1):n*o, j)'.*x + q(j));
          y = sum(x.*p(st, :)'+repmat(q', [1 size(x, 2)]), 2);
%           y = softmax(y);
%           u(i+n*(o-1),j) = exp(-0.5*(((x(i-n*(o-1))-c(i,j))/s(i,j))^2));
%           u(1+n*(o-1):n*o,j) = exp(-0.5.*(((x-c(1+n*(o-1):n*o,j)')./s(1+n*(o-1):n*o, j)').^2));
          u = exp(-0.5.*(((x'-c(st,:))./s(st, :)).^2));

%           w(j) = w(j) * u(i-n*(o-1),j);
%            w(j) = prod(u(1+n*(o-1):n*o,j));
          w = prod(u);
%           w = w./sum(w);
%        end
%        a = a + w(j)*y(j);
       a = w*y;
%        b = b + w(j);
       b = sum(w);
%     end
    ys = a/b;
%     ys = a;
end
