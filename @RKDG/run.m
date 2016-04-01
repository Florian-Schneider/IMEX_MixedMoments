function run( obj )
%RUN Wrapper for a complete run
if usejava('jvm') && feature('ShowFigureWindows')
    h = waitbar(0,'Calculating');
else
    msg = sprintf('t = %2.3f, Estimated time: %10.2f seconds',obj.t(1),0);
    fprintf(msg)
end

tic0 = tic;
dt = obj.dt;

last_save = tic0;

A = 0; b = 1; c = 0; %Butcher tableau for Euler step


obj.MomentModel.TimeStepping.A = A;
obj.MomentModel.TimeStepping.b = b;
obj.MomentModel.TimeStepping.dt = dt;
for n=obj.t_cnt+1:size(obj.t,2)-1
    obj.t_cnt = n;
    
    if ismember(obj.t_cnt,obj.t_Frame_ind)
        obj.u_Frame(:,:,:,obj.Frame_cnt) = obj.u;
        runtime = toc(tic0);
        estimated_time = round((size(obj.t,2)-n)*runtime/n);
        if usejava('jvm') && feature('ShowFigureWindows')
            waitbar(n/(size(obj.t,2)-1),h,sprintf('Calculating - Time left: ~%d seconds\n',estimated_time));
        else
            fprintf(repmat('\b',1,length(msg)));
            msg = sprintf('t = %2.3f, Estimated time: %10.2f seconds',obj.t(n)+dt,estimated_time);
            fprintf(msg)
        end
        obj.Frame_cnt = obj.Frame_cnt +1;
    end
    
    k = zeros([size(obj.u),length(A)]);
    for i=1:length(A) %RK
        obj.MomentModel.TimeStepping.Stage = i;
        if i==1
            ustage = obj.u;
        else
            ustage = obj.u+dt*sum(bsxfun(@times,permute(A(i,:),[1 3 4 2]),k),4);
        end
        uB = ustage;
        k(:,:,:,i) = DifferentialOperator(obj,uB,obj.t(n)+c(i)*dt);
        if n==1 && i==1
            obj.u0_Frame(:,:,1) = obj.u_0_rec;
        end
    end
    ustage = obj.u+dt*sum(bsxfun(@times,permute(b,[2 3 4 1]),k),4); %Update fluxes
    ustage = DifferentialOperatorRHS(obj,ustage,obj.t(n)); %Implicit part
    obj.u = ustage;
    
    if obj.save_flag && toc(last_save)>save_table(toc(tic0)/n)
        obj.save(obj.Path,obj.Filename);
        last_save = tic;
        disp('saved')
        if toc(tic0)>obj.destroytime
            eval(obj.restartCode);
            fprintf('I killed myself after %f hours of running\n',toc(tic0)/3600);
            exit;
        end
    end
    
end




obj.runtime = toc(tic0);
obj.t_cnt= obj.t_cnt+1;
if ismember(obj.t_cnt,obj.t_Frame_ind)    
    obj.u_Frame(:,:,:,obj.Frame_cnt) = obj.u;
end

if usejava('jvm') && feature('ShowFigureWindows')
    delete(h);
else
    fprintf(repmat('\b',1,length(msg)));
    disp('Done');
end

if obj.save_flag
    obj.save(obj.Path,obj.Filename);
    disp('saved')
end

end


function t = save_table(s)
t = 3600;
if s<=120
    t = 1200;
end
if s<=60
    t = 600;
end
if s<=10
    t = 60;
end

end

