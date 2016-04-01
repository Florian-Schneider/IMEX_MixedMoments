function printCuts(R,filename,x,y,t)

DichteCut = R.evaluateDensity(x,y,t);
mkdir(filename)
fid = fopen([filename class(R.MomentModel),'_' num2str(R.MomentModel.MomentOrder) '-Q' num2str(R.MomentModel.QuadratureOrder) '.txt'],'w+');
fprintf(fid,'x y');
fprintf(fid,' t%d',1:length(t));
fprintf(fid,'\n');

fprintf(fid,[repmat('%f ',1,length(t)+2) '\n'],[x,y,DichteCut]');
fclose(fid);
end