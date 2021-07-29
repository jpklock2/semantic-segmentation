function [ys, cent, sig, pi, qi] = CANFIS(xt, ydt, nEp, alfa, plotFigs, inFis)

% close all;
% clear all;
% clc;

% load ydt4;
% load ydv4;

% ydt = ydt4;
% ydv = ydv4;
% xt = xt4;
% xv = xv4;

% [~, rawresult] = system('node C:\Users\rodri\OneDrive\Documentos\UFMG\PFC\MachineLearning_v4\index.js');
% resNode = jsondecode(rawresult);
% ydv2 = resNode.y;

% resNode.all.x(resNode.all.x <= 1e-3) = 0;

% ydt = classesOri;
% ydv = ;
% xt = pixelsOri;
% xv = resNode.all.x;

nRe = size(ydt, 2);
nO = size(ydt, 2);
% alfa=0.0001;
nMf=length(xt(1,:));
% nEp=100;
npt=length(ydt);
% npv=length(ydv);

xmin = min(xt);
xmax = max(xt);
delta=(xmax-xmin)/(nRe);

cent = zeros(nMf, nRe);
sig = zeros(nMf, nRe);
pi = zeros(nMf, nRe, nO);
qi = zeros(1, nRe);
% cnt = 0;

% for o = 1:nO
%     cnt = cnt+b*nRe;
    for j=1:nRe
       for i=1:nMf
           if nargin > 4
                cent(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(2);
                sig(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(1);
           else
                cent(i,j) = xmin(i) + (j-1)*delta(i) + delta(i)/2;
                sig(i,j) = 0.1+rand; %0.5*delta(i)*sqrt(1/(2*log(2)));
           end
          pi(i,j,:) = rand([1 1 nO]);
       end
       qi(j) = rand;
    end
% end

fprintf('\n');
fP = 0;
cost = zeros(nEp, 1);

for l=1:nEp
    
    if l > 1 && mod(l, 50) == 0 && l < nEp
        alfa = alfa/3;
    end
    
    prop = max(round(25*(l/nEp)), 1);
    if prop >= fP
        if l > 1
            fprintf(repmat('\b', 1, length(str1) + length(str2) - 3));
        end
        str1 = ['Epoch = ' num2str(l) '/' num2str(nEp) ' (' num2str(round(100*(l/nEp),2)) '%%)\n'];
        fprintf(str1);
        str2 = ['[' repmat(char(9608), 1, prop) repmat(char(8195), 1, 25-prop) ']\n'];
        fprintf(str2);
        fP = min(fP + 25*(1/nEp), 25);
    end
    
    ys = zeros(npt, nO);
    for k=1:npt
        y = zeros(nO, nO);
        w = zeros(nO, nO);
        b = zeros(nO, 1);
        for o=1:nO
            %calcular ys, y, w e b
            [ys(k,o),y(:,o),w(o,:),b(o)] = saida(xt(k,:),pi(:,:,o),qi,sig,cent,nRe,nMf,1);
        end
%         ys(k,:) = softmax(ys(k,:)')';
        if any(isnan(ys(k,:)))
            disp(k);
            dbgg = 1;
        end
        
        st = 1:nMf;
        dedys = ys(k,:) - ydt(k,:);
        dysdwj = (y - ones(nO,1)*ys(k,:)) / b(1);
%         dedys = sum(ys(k,:) - ydt(k,:));
%         dysdwj = sum((y - ones(nO,1)*ys(k,:)) / b(1), 2);
        dysdyj = w(1,:) ./ b(1);
        dwdcij = w(1,:).*((xt(k,:)'-cent(st,:))./(sig(st,:).^2));
        dwdsij = w(1,:).*(((xt(k,:)'-cent(st,:)).^2)./(sig(st,:).^3));
        dydqj = 1;
        dedcij = dwdcij.*(sum(dedys'.*dysdwj, 2))';
        dedsij = dwdsij.*(sum(dedys'.*dysdwj, 2))';
%         dedcij = dwdcij.*(dedys.*dysdwj)';
%         dedsij = dwdsij.*(dedys.*dysdwj)';
        dedqij = dedys.*dysdyj*dydqj;
        cent(st,:) = cent(st,:) - (alfa.*dedcij);
        sig(st,:) = sig(st,:) - (alfa.*dedsij);
        qi = qi - (alfa*nMf*dedqij);
        
        for o=1:nO
%             st = 1:nMf;
            dedys = ys(k,o) - ydt(k,o);
%             for j=1:nRe
%                 dysdwj = (y(j) -  ys(k,o)) / b;
%                 dysdwj = (y(:,o) -  ys(k,o)) / b(o);
%                 dysdyj = w(j) / b;
                dysdyj = w(o,:) ./ b(o);

%                 for i=1+nMf*(o-1):nMf*o
%                     dwdcij = w(j)*((xt(k,i-nMf*(o-1))-cent(i,j))/(sig(i,j)^2));
%                     dwdcij = w(j).*((xt(k,:)'-cent(st,j))./(sig(st,j).^2));
%                     dwdcij = w(o,:).*((xt(k,:)'-cent(st,:))./(sig(st,:).^2));
%                     dwdsij = w(j)*(((xt(k,i-nMf*(o-1))-cent(i,j))^2)/(sig(i,j)^3));
%                     dwdsij = w(j).*(((xt(k,:)'-cent(st,j)).^2)./(sig(st,j).^3));
%                     dwdsij = w(o,:).*(((xt(k,:)'-cent(st,:)).^2)./(sig(st,:).^3));
%                     dydpij = xt(k,i-nMf*(o-1));
                    dydpij = xt(k,:)';
%                     dedcij = dedys*dysdwj.*dwdcij;
%                     dedcij = dwdcij.*(dedys.*dysdwj)';
%                     dedsij = dedys*dysdwj.*dwdsij;
%                     dedsij = dwdsij.*(dedys*dysdwj)';
%                     dedpij = dedys*dysdyj.*dydpij;
                    dedpij = (dedys.*dysdyj).*dydpij;

%                     cent(st,:) = cent(st,:) - (alfa.*dedcij);
%                     sig(st,:) = sig(st,:) - (alfa.*dedsij);
                    pi(:,:,o) = pi(:,:,o) - (alfa.*dedpij);
%                 end
%             end
        end
    end
    
    cost(l) = mean(sum((1/2).*((ys - ydt).^2))./length(ys));
%     cost(l) = mean(abs(ys-ydt));
    
    if (l > 5 && plotFigs)
        children = get(gca, 'children');
        if (~isempty(children)); delete(children(1)); end
        plot(5:nEp, cost(5:end), 'r*');
        title('ANFIS Training Cost');
        xlabel('Nยบ Epochs');
        ylabel('Cost');
        pause(0.1);
    end
    
end

fprintf('\n');

% for k=1:npt
%     for o=1:nO
%         [YsOut(k, o),y,w,b] = saida(pixelsOri(k,:),pi,qi,sig,cent,nRe,nMf,o);
%     end
% end

% for k=1:npv
%     [ysaida(k),y,w,b] = saida(xv(k,:),p,q,s,c,nRe,nMf);
% end
% subplot(131); 

% ysaida = (ysaida - min(ysaida))./(max(ysaida) - min(ysaida));
% plot(ydv); hold on; plot(ysaida);
% legend({'Original', 'MATLAB'});
% figure
% %[~, rawresult] = system('node C:\Users\joao.ferreira\Desktop\MachineLearning_v2\index.js');
% %ydv2 = jsondecode(rawresult);
% ydv2 = (ydv2 - min(ydv2))./(max(ydv2) - min(ydv2));
% plot(ydv); hold on; plot(ydv2);
% legend({'Original', 'Node'});
end