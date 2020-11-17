function framedata = emma_2d(folderstr,type,nodenum,index,indx,indy,ftype)

% function framedata = emma_2d(folderstr,type,nodenum,index,indx,indy)

nodeindex = nextpow2(nodenum);

if (strcmp(type,'edensity')==1)
    typestr = '/DENSITY/electron/echarge-';
elseif (strcmp(type,'idensity')==1)
    typestr = '/DENSITY/ion/icharge-';
elseif (strcmp(type,'e1')==1)
    typestr = '/FLD/e1/e1-';
elseif (strcmp(type,'e2')==1)
    typestr = '/FLD/e2/e2-';
elseif (strcmp(type,'e3')==1)
    typestr = '/FLD/e3/e3-';
elseif (strcmp(type,'b1')==1)
    typestr = '/FLD/b1/b1-';
elseif (strcmp(type,'b2')==1)
    typestr = '/FLD/b2/b2-';
elseif (strcmp(type,'b3')==1)
    typestr = '/FLD/b3/b3-';
elseif (strcmp(type,'ej1')==1)
    typestr = '/CURRENT/electron/j1/j1-';
elseif (strcmp(type,'ej2')==1)
    typestr = '/CURRENT/electron/j2/j2-';
elseif (strcmp(type,'ej3')==1)
    typestr = '/CURRENT/electron/j3/j3-';
elseif (strcmp(type,'tj1')==1)
    typestr = '/CURRENT/total/j1/j1-';
elseif (strcmp(type,'tj2')==1)
    typestr = '/CURRENT/total/j2/j2-';
elseif (strcmp(type,'tj3')==1)
    typestr = '/CURRENT/total/j3/j3-';
else
    disp('type not recognized.');
    disp('type = edensity or idensity or e1 or e2 or e3');
    
end

file2str = '.dat';
framedata = zeros(2^indx,2^indy);

for node = 0 : 2^nodeindex - 1
        filestr = [folderstr typestr getnumstr(index,8) '-'...
                            getnumstr(node,4) file2str];
        if (ftype==1)
            fid = fopen(filestr);
            framein = fread(fid,'double');
            frame = reshape(framein, 2^indx,2^(indy-nodeindex));
            framedata(:, node*2^(indy-nodeindex)+1:(node+1)*2^(indy-nodeindex)) = ...
                frame;
            fclose(fid);
        elseif (ftype==2)
            framein=textread(filestr);
            frame = reshape(framein, 2^indx,2^(indy-nodeindex));
            framedata(:, node*2^(indy-nodeindex)+1:(node+1)*2^(indy-nodeindex)) = ...
                frame;
        end
end
