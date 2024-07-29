function[xp,yp] = linehist(x,y,varargin)
% linehist(x,y,varargin)

% INPUTS
%   x: bin centers (evenly spaced)
%   y: data
%
% (Optional)
%   'normalize': double - normalizes y axis to given value. E.g.
%                ['normalize', 100] will make all histograms add to 100%
%   'Color': char or 3-element vector - specifies plot color
%   'LineType': char - Changes line plot type (e.g. ':' or '--'); 
%   'suppress_plot': boolean - suppresses plotting when true
%   'fill': double between [0 1] - Specifies the alpha value of the fill. 

normcount = 0;
color2plot = [];
make_plot = true; 
linetype = '-'; 
fillplot = 0; 
for i = 1:length(varargin)
    if strcmp(varargin{i},'normalize')
        normcount = varargin{i+1};
    elseif ischar(varargin{i}) && length(varargin{i})==1 && ~strcmp(varargin{i-1},'LineType')
        color2plot = varargin{i};
    elseif strcmp(varargin{i},'Color')
        color2plot = varargin{i+1};
    elseif strcmp(varargin{i},'LineType')
        linetype = varargin{i+1};
    elseif strcmp(varargin{i},'suppress_plot')
        make_plot = false; 
    elseif strcmp(varargin{i},'fill')
        fillplot = varargin{i+1}; 
    end
end

if numel(y)>1
    x = reshape(x,1,[]);
    dx = diff(x(1:2));
    xedge = [x(1)-.5*dx x+.5*dx];
    xedge2 = [-inf xedge inf]; 
    
    for i = 1:size(y,2)
        n = histcounts(y(:,i),xedge2); 
        if ( n(1)+n(end) ) > 0
            not_inclusive = true; 
            warning('Histogram edges do not contain all data'); 
        else
            not_inclusive = false; 
        end
        n(2) = n(1)+n(2); n(end-1) = n(end-1)+n(end); n = n(2:(end-1)); 
        b = x;

        d = mode(diff(b));
        drep = [-(d/2)*ones(1,length(b)); (d/2)*ones(1,length(b))];

        xp = reshape(repmat(reshape(b,1,[]),2,1)+drep,1,[]) ;
        yp = [0 reshape(repmat(reshape(n,1,[]),2,1),1,[]) 0];

        if normcount > 0
            yp = yp./sum(n)*normcount;
        end

        xp = [xp(1) xp(1:end) xp(end)];

        if make_plot
            hold on; 
            if ~isempty(color2plot)
                ph = plot(xp,yp,linetype,'Color',color2plot,'LineWidth',2);
            else
                ph = plot(xp,yp,linetype,'LineWidth',2);
            end
            if not_inclusive
                yl = ylim; 
                text(xp(end),yl(1),'+','FontSize',16,'VerticalAlignment','top'); 
                text(xp(1),yl(1),'-','FontSize',16,'VerticalAlignment','top'); 
            end
    
            if fillplot>0
                patch(ph.XData,ph.YData,ph.Color,'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off'); 
            end

            xlim([xp(1) xp(end)]);
        end
    end
else
    
    [n,b] = hist(x,y);
    n = n'; 
    d = mode(diff(b));
    drep = [-(d/2)*ones(1,length(b)); (d/2)*ones(1,length(b))];
    
    xp = reshape(repmat(reshape(b,1,[]),2,1)+drep,1,[]) ;
    xp = [xp(1) xp(1:end) xp(end)];
    for j = 1:size(n,2)
        yp = [0 reshape(repmat(reshape(n(:,j),1,[]),2,1),1,[]) 0];
        if normcount
            yp = yp./sum(yp);
        end

        if make_plot
            hold on; 
            if ~isempty(color2plot)
                ph = plot(xp,yp,linetype,'Color',color2plot,'LineWidth',2);
            else
                ph = plot(xp,yp,linetype,'LineWidth',2);
            end

            if fillplot>0
                patch(ph.XData,ph.YData,ph.Color,'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off'); 
            end

            xlim([xp(1) xp(end)]);
        end
    end
end
    