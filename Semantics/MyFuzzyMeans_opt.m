function [idx2, U, J, centroids] = MyFuzzyMeans_opt(X, K)

%------------------------------------------------------------------------
% Algoritmo Fuzzy C-Means Otimizado

plots = 0; % plot em 2 dimensões
plots3 = 0; % plot em 3 dimensões
prints = 0; % print dos valores dos resultados

% definicao das variaveis numericas
m = 2; % numero M da matriz U
iter = 1; % iteracao inicial
Nmax = 101; % numero maximo de iteracoes permitido
tol = 0.0001; % tolerancia para parar estibilizacao da funcao objetivo

%% criacao da matriz aleatoria U
n = size(X,1); % numero de elementos de X
rnd = rand(n, K); % cria uma matriz de elementos aleatorios
U = rnd./(sum(rnd, 2)); % normaliza matriz por linhas
U_sq = U.^2; % calcula matriz U quadrada para diminuir calculos

%% achando os centroides iniciais
centroids = ((X'*U_sq)./sum(U_sq))'; % calculo dos centroides por cluster

%% encontrando valor inicial da funcao objetivo
dist_matrix = zeros(n, K); % cria matriz de distancias
for j = 1:K % iterando pelos centroides (clusters)
    dist_matrix(:, j) = sqrt(sum(((X-ones(n, 1).*centroids(j, :)).^2), 2)); % calcula a distancia de todos os elementos ao cluster
end
J(iter) = sum(sum((dist_matrix.^2).*U_sq)); % calcula a funcao objetivo

%% caso se deseje ver os valores da funcao objetivo
if prints
    fprintf('\nFunção objetivo inicial: %.5f\n', J(iter));
end

% abaixo, caso se deseje plotar os clusters
if plots
    figure(1);
    hold on;
    plot(X(:,1), X(:,2),'b.');
    grid on;
    xdata = centroids(:,1);
    ydata = centroids(:,2);
    h = plot(xdata, ydata, 'ko', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize', 10);
    pause(1)
end

% abaixo, caso se deseje plotar os clusters em 3 dimensoes
if plots3
    figure(1); subplot(1, 2, 2);
    title('Pixels in the RGB space')
    plot3(X(:,1),X(:,2),X(:,3),'b.');
    axis([0 1 0 1 0 1]);
    hold on;
    axis square; grid on;
    xdata = centroids(:,1);
    ydata = centroids(:,2);
    zdata = centroids(:,3);
    h = plot3(xdata, ydata, zdata, 'ko', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize', 10);
    pause(1)
end

%% loop principal
change = Inf; % valor inicial da deteccao de mudancas
% itera ate o numero maximo de iteracoes ou a funcao objetivo mudar menos que a tolerancia
while (change > tol && iter < Nmax)
    
    % calcula os centroides
    centroids = ((X'*U_sq)./sum(U_sq))'; % calculo dos centroides por cluster
    
    % plotando novos centroides em 2 dimensoes
    if plots
        xdata = centroids(:,1);
        ydata = centroids(:,2);
        pause(1)
        set(h,'XData',xdata,'YData',ydata);
    end
    
    % plotando novos centroides em 3 dimensoes
    if plots3
        xdata = centroids(:,1);
        ydata = centroids(:,2);
        zdata = centroids(:,3);
        pause(1)
        set(h,'XData',xdata,'YData',ydata,'ZData',zdata);
    end
    
    % calculo da nova U
    dist_matrix = zeros(n, K); % cria matriz de distancias
    for j = 1:K % iterando pelos centroides (clusters)
        dist_matrix(:, j) = sqrt(sum(((X-ones(n, 1).*centroids(j, :)).^2), 2)); % calcula a distancia de todos os elementos ao cluster
    end
    dist_sq = dist_matrix.^(-2/(m-1)); % calcula a matriz de distancias ao quadrado
    U = dist_sq./sum(dist_sq, 2); % calcula novo U
    U_sq = U.^2; % calcula matriz U quadrada
        
	% calculo da funcao objetivo
    iter = iter + 1; % incrementa a iteracao
    J(iter) = sum(sum((dist_matrix.^2).*U_sq)); % calcula a funcao objetivo
    
    % verificacao do criterio de parada
    if prints % caso se deseje printar o criterio de parada por iteracao
        fprintf('Função objetivo iteração %d: %.5f\n', iter-1, J(iter));
    end
    change = abs(J(iter)-J(iter-1)); % calcula novo valor de mudanca
    
end

[~, idx2] = sort(U, 2);
idx2 = idx2(:, end);

% plotando o resultado final do Fuzzy C-means para 2 dimensoes
if plots
    clus = 1:K; % define clusters
    colors = {'b.', 'r.', 'c.', 'm.', 'y.', 'k.'};  % define cores padrao
    for i = 1:K % itera pelos clusters
        [~, idx] = sort(U, 2); % ordena U por linha
        indexes = find(idx(:, end)==clus(i)); % pega o maior indice de cada linha ordenada (cluster ao qual o elemento pertence)
        plot(X(indexes,1), X(indexes,2), colors{i}); % plota cluster na determinada cor
    end
    h = plot(xdata, ydata, 'ko', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize', 10); % plota clusters

    % plotando a funcao objetivo
    figure(2);
    hold on;
    plot(1:iter, J(1:iter), 'b--', 'LineWidth', 2);
    grid on;
end

% plotando o resultado final do Fuzzy C-means para 3 dimensoes
if plots3
    clus = 1:K; % define clusters
    colors = {'b.', 'r.', 'c.', 'm.', 'y.', 'k.'};  % define cores padrao
    for i = 1:K % itera pelos clusters
        [~, idx] = sort(U, 2); % ordena U por linha
        indexes = find(idx(:, end)==clus(i)); % pega o maior indice de cada linha ordenada (cluster ao qual o elemento pertence)
        plot3(X(indexes,1), X(indexes,2), X(indexes,3), colors{i}); % plota cluster na determinada cor
    end
    h = plot3(xdata, ydata, zdata, 'ko', 'LineWidth', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor','g', 'MarkerSize', 10); % plota clusters

    % plotando a funcao objetivo
    figure(2);
    hold on;
    plot(1:iter, J(1:iter), 'b--', 'LineWidth', 2);
    grid on;
end
