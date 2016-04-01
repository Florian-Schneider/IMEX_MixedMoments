function plot( obj , varargin)
%plot plot-function for a Kinetic Scheme

p = inputParser;
p.addParameter('Scaling','lin');
p.addParameter('caxis',[]);
p.addParameter('axis',[]);
p.addParameter('Time',[-inf,inf]);
p.addParameter('Frames',inf);
p.addParameter('colorbar',true);
p.addParameter('colormap',parula(1000));
p.addParameter('title',true);
p.addParameter('xticks',[]);
p.addParameter('xlabel','$x$');
p.addParameter('yticks',[]);
p.addParameter('ylabel','$y$');
p.addParameter('handle',[]);
p.addParameter('flat',true);
p.KeepUnmatched = true;
p.parse(varargin{:});

if ~isempty(fields(p.Unmatched))
    F = fields(p.Unmatched);
    for i=1:length(fields(p.Unmatched))
        warning(sprintf('It seems that there is an unmatched parameter-value pair: %s',F{i})); %#ok<*SPWRN>
    end
end
Time = [max(min(p.Results.Time),min(obj.t_Frame)),min(max(p.Results.Time),max(obj.t_Frame))];
Frames = min(p.Results.Frames,obj.Frames);

u0 = permute(obj.MomentModel.Density(permute(obj.u_Frame,[3 1 2 4])),[2 3 4 1]);
%u0 = permute(obj.u_Frame(:,:,1,:),[1 2 4 3]);
u0 = u0/max(obj.MomentModel.BasisFunctionsAtQuadratureSqueezed(1,:));
scale_min = min(reshape(u0(:,1,:),[],1));
scale_max = max(reshape(u0(:,1,:),[],1));

if isempty(p.Results.caxis)
    CA = [scale_min,scale_max];
else
    CA = p.Results.caxis;
end



if isempty(p.Results.axis)
    A = [obj.Grid.xmin,obj.Grid.xmax,obj.Grid.ymin,obj.Grid.ymax];
else
    A = p.Results.axis;
end

CB = p.Results.colorbar;
scaling = p.Results.Scaling;

if isempty(p.Results.handle)
    handle = figure(27);
    clf;
    set(handle,'Renderer','painters')
else
    handle = p.Results.handle;
    subplot(handle);
    axes(handle);
    cla;
end

I = unique(interp1(obj.t_Frame,1:obj.Frames,linspace(Time(1),Time(2),Frames),'nearestneighbor'));

for i=I
    if strcmp(scaling,'lin') || strcmp(scaling,'nlin')
        logarithmic = false;
        
        if strcmp(scaling,'lin')            
            obj.Grid.plot(min(max(u0(:,:,i),CA(1)),CA(2)),p.Results.flat,logarithmic);            
            caxis(CA)
        else
            obj.Grid.plot(u0(:,:,i),p.Results.flat,logarithmic);
        end
    else
        logarithmic = true;
        obj.Grid.plot(u0(:,:,i),p.Results.flat,logarithmic);
        if strcmp(scaling,'log')
            caxis(log10([max(CA(1),1e-10),CA(2)]))
        end
    end
    xlabel(p.Results.xlabel,'Interpreter','latex','FontSize',13);
    ylabel(p.Results.ylabel,'Interpreter','latex','FontSize',13);
    if islogical(p.Results.title) && p.Results.title
        title(sprintf('t = %f\n',obj.t_Frame(i)));
    elseif ischar(p.Results.title)
        title(p.Results.title)
    end

    shading interp
    if CB
        colorbar
    end
    axis(A);
    if p.Results.flat
        view(2)
    else
        view(3)
    end
    if ~isempty(p.Results.xticks)
        set(gca,'XTick',p.Results.xticks);
    end
    if ~isempty(p.Results.yticks)
        set(gca,'YTick',p.Results.yticks);
    end
    
    drawnow;
end


end