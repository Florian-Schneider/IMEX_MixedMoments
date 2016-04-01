function printSolutions( K,varargin)

p = inputParser;
p.addParamValue('folder','Images/Dummy');
p.addParamValue('caxis',[]);
p.addParamValue('axis',[]);
p.addParamValue('Time',[]);
p.addParamValue('Scaling','nlin');
p.addParamValue('colormap',parula(1000));
p.addParamValue('xticks',[]);
p.addParamValue('xlabel','$x$');
p.addParamValue('yticks',[]);
p.addParamValue('ylabel','$y$');
p.addParamValue('modelname',[]);

p.KeepUnmatched = true;
p.parse(varargin{:});
if ~isempty(fields(p.Unmatched))
    F = fields(p.Unmatched);
    for i=1:length(fields(p.Unmatched))
        warning(sprintf('It seems that there is an unmatched parameter-value pair: %s',F{i})); %#ok<*SPWRN>
    end
end

if ~iscell(K)
    tmp{1} = K;
    K = tmp; %put K into a 1x1 cell array if it is only a single solution file
end
folder = p.Results.folder;
if ~isdir(folder)
    mkdir(folder);
end


if isempty(p.Results.caxis)
    scale_min = inf;
    scale_max = -inf;
else
    scale_min = p.Results.caxis(:,1);
    scale_max = p.Results.caxis(:,2);
end
if isempty(p.Results.caxis)
    for j=1:length(p.Results.Time)
        for i=1:length(K)
            [~,I] = min(abs(K{i}.t_Frame-p.Results.Time(j)));
            u0 = K{i}.u0_Frame(:,:,I);
            u0 = u0/max(K{i}.MomentModel.BasisFunctionsAtQuadratureSqueezed(1,:));
            scale_min = min(min(u0(:)),scale_min);
            scale_max = max(max(u0(:)),scale_max);
        end
    end
end
for j=1:length(p.Results.Time)
    for i=1:length(K)
        Param.Scaling = p.Results.Scaling;
        if size(scale_min,1)==length(p.Results.Time)
            Param.caxis = [scale_min(j),scale_max(j)];
        else
            Param.caxis = [scale_min,scale_max];
        end
        Param.axis = p.Results.axis;
        
        Param.Time = p.Results.Time(j)*[1,1];
        Param.colorbar = false;
        Param.colormap = p.Results.colormap;
        Param.title = false;
        Param.xticks = p.Results.xticks;
        Param.xlabel = p.Results.xlabel;
        Param.yticks = p.Results.yticks;
        Param.ylabel = p.Results.ylabel;
        K{i}.plot(Param);
        filename = [folder '/' class(K{i}.MomentModel),'_' num2str(K{i}.MomentModel.MomentOrder),'-' p.Results.modelname '-Frame' num2str(j,'%04d')];
        set(gcf,'Color','w');
        generate2Dplot(gcf,filename,'colorbar.tex')
        %export_fig(filename,'-pdf','-m4','-ZBuffer')
    end
end

% fid = fopen([folder '/time.txt'],'w+');
% fprintf(fid,'%f\n',p.Results.Time);
% fclose(fid);
% 
% fid = fopen([folder '/colorbar.tex'],'w+');
% switch p.Results.Scaling
%     case 'nlin'
%         fprintf(fid,'\\begin{tikzpicture}\n \\begin{axis}[\n hide axis,\n scale only axis,\n height=0pt,\n width=0pt,\n colormap/jet,\n colorbar horizontal,\n point meta min=%f,\n point meta max=%f,\n colorbar style={\n width=10cm}]\n\\addplot [draw=none] coordinates {(0,0)};\n\\end{axis}\n\\end{tikzpicture}',scale_min,scale_max);
%     case 'nlog'
%         fprintf(fid,'\\begin{tikzpicture}\n \\begin{axis}[\n hide axis,\n scale only axis,\n height=0pt,\n width=0pt,\n colormap/jet,\n colorbar horizontal,\n point meta min=%f,\n point meta max=%f,\n colorbar style={\n width=10cm}]\n\\addplot [draw=none] coordinates {(0,0)};\n\\end{axis}\n\\end{tikzpicture}',log10(scale_min),log10(scale_max));
% end
% %fprintf(fid,'\\begin{tikzpicture}\n \\begin{axis}[\n hide axis,\n scale only axis,\n height=0pt,\n width=0pt,\n colormap/jet,\n colorbar horizontal,\n point meta min=%f,\n point meta max=%f,\n colorbar style={\n width=10cm,\n xtick={18,20,25,...,45}\n}]\n\\addplot [draw=none] coordinates {(0,0)};\n\\end{axis}\n\\end{tikzpicture}',scale_min,scale_max);
% fclose(fid);

%class(K.Case);
end

