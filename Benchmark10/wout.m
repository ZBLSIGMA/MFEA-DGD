 load result;

for i = 1 : 10
    filename=['MTOCSO_P' num2str(i) '.txt'];
    fid=fopen(filename,'w');
    for r = 10 : 10 : 1000
        fprintf(fid,'%d,',r*200,data_MFEA(i).EvBestFitness(1:59,r));
        fprintf(fid,'%d',data_MFEA(i).EvBestFitness(60,r));
        fprintf(fid,'\r\n');
    end
    fclose(fid);
end

