
snow=[245/255,254/255,253/255]; % grama clara
yellow=[253/255,230/255,75/255]; % grama normal
orange=[237/255,112/255,20/255]; % grama escura
red=[208/255,49/255,45/255]; % floresta clara
pink=[246/255,154/255,205/255]; % floresta escura
purple=[163/255,44/255,196/255]; % cimento (?)
blue=[58/255,67/255,186/255]; % terra cinza
green=[59/255,177/255,67/255]; % terra marrom
grey=[108/255,98/255,109/255]; % pouca mata
brown=[35/255,23/255,9/255]; % asfalto
chartreuse=[172/255,252/255,56/255]; % estradas
tortilla=[154/255,123/255,79/255]; % mata cinza
mint=[152/255,237/255,195/255]; % rio
hotpink=[255/255,22/255,149/255]; % pedra
% silver=[173/255,173/255,199/255]; 
% artic=[130/255,237/255,253/255];

myColors = [snow; yellow; orange; red; pink; ...
            purple; blue; green; brown; grey; ...
            chartreuse; tortilla; mint; hotpink]; %; ...
            % silver; artic];
        
myColorsNames = {'snow'; 'yellow'; 'orange'; 'red'; 'pink'; ...
                'purple'; 'blue'; 'green'; 'brown'; 'grey'; ...
                'chartreuse'; 'tortilla'; 'mint'; 'hotpink'}; %; ...
                % 'silver'; 'artic'};

% myColorsClasses = {'Light Grass'; 'Common Grass'; 'Dark Grass'; 'Light Forest'; 'Dark Forest'; ...
%                 'Gray Soil'; 'Grey Earth'; 'Brown Earth'; 'Little Foliage'; 'Asphalt Road'; ...
%                 'Earth Road'; 'Grey Forest'; 'River'; 'Stone'}; %; ...
%                 
% fig = figure; hold on;
% x1 = [1:5 5*ones(1, 3) 5:-1:1 ones(1, 3)];
% y1 = [6*ones(1, 5) 7:9 10*ones(1, 5) 9:-1:7];
% xpos = 3;
% ypos = 10.5;
% for cor = 1:14
%     fill(x1, y1, myColors(cor, :));
%     x1 = x1+5;
%     text(xpos, ypos, [num2str(cor) ' - ' myColorsClasses{cor}], 'FontSize', 18, 'HorizontalAlignment', 'center');
%     xpos = xpos+5;
%     if cor == 7
%         x1 = [1:5 5*ones(1, 3) 5:-1:1 ones(1, 3)];
%         y1 = y1-5;
%         xpos = 3;
%         ypos = ypos-5;
%     end
% end
% axis([0 36 0 11]);
% 
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% ax = gca;
% ax.Visible = 'off';
% 
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'D:\TCC\TCC\figs\colors.pdf','-dpdf','-r0');
% 
% dbg = 1;
