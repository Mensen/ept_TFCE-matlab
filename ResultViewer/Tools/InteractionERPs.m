AvgF11 = squeeze(mean(cell2mat(Data(1,1)')));
AvgF12 = squeeze(mean(cell2mat(Data(2,1)')));
AvgF21 = squeeze(mean(cell2mat(Data(1,2)')));
AvgF22 = squeeze(mean(cell2mat(Data(2,2)')));

IdCh = 55;

F11  =  AvgF11(IdCh,:);
F12  =  AvgF12(IdCh,:);
F21  =  AvgF21(IdCh,:);
F22  =  AvgF22(IdCh,:);

plot(F11,'DisplayName','F11','YDataSource','F11');hold all;plot(F12,'DisplayName','F12','YDataSource','F12');plot(F21,'DisplayName','F21','YDataSource','F21');plot(F22,'DisplayName','F22','YDataSource','F22');hold off;figure(gcf);


% Topoplot at time

x=Data(:);
x=squeeze(mean(cell2mat(x)));

s = 88;

x1 = x(:,s);
figure
topoplot(x1,e_loc...
    ,'electrodes', 'on'... %display markers ("labels" shows electrode names
    ,'whitebk', 'on' ...
    ,'colormap', colormap(jet(64)) ...
    ,'shading', 'interp'... % (flat or interp) useless with style is "fill"
    ,'style', 'both'...
    ,'plotrad', max(abs(cell2mat({e_loc.radius})))); % 'map' 'contour' 'both' 'fill' 'blank'
