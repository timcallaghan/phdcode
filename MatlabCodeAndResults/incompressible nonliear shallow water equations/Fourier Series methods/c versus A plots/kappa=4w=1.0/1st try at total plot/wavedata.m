clear all
load AmpsandC


T = 30./AmpsandC(1,:);

max = length(AmpsandC);

fid = fopen('wavedata.txt','wt');
fclose(fid);
fid = fopen('wavedata.txt','at');

for i = 1:max
    fprintf(fid,'%.16e',AmpsandC(1,i));
    fprintf(fid,' ');
end

for i = 1:max
    fprintf(fid,'%.16e',AmpsandC(2,i));
    fprintf(fid,' ');
end

for i = 1:max
    fprintf(fid,'%.16e',AmpsandC(3,i));
    fprintf(fid,' ');
end

for i = 1:max
    fprintf(fid,'%.16e',AmpsandC(4,i));
    fprintf(fid,' ');
end

for i = 1:max-1
    fprintf(fid,'%.16e',T(1,i));
    fprintf(fid,' ');
end
fprintf(fid,'%.16e',T(1,max));

fclose(fid)