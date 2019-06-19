%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Library load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DB      = importdata('materialDB.csv'          ,',');
DB_epm  = importdata('materialDB_epm.csv'      ,',');
DB_TB_V = importdata('materialDB_TBVogl.csv'   ,',');
DB_TB_K = importdata('materialDB_TBKlimeck.csv',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M{1}=Material;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Here is a patch to be able to load the tables in Matlab AND Octave %%%%%%
% Matlab see the header in multiple cells while Octave see the header in one cell only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB.textdata(1,:))==1
  
    DB.textdata{1,1}=[DB.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB.textdata{1,1},',');
    idx=strfind(DB.textdata{1,1},[',' M{1} ',']);
    idxM=find(idxM==idx);
    
    M{2} = DB.data(:,idxM);
else

    for i=1:length(DB.textdata(1,:))
      idx=strcmp(DB.textdata{1,i},M{1});
      if idx==1
        M{2}=DB.data(:,i-1);
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB_epm.textdata(1,:))==1
    
    DB_epm.textdata{1,1}=[DB_epm.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB_epm.textdata{1,1},',');
    idx=strfind(DB_epm.textdata{1,1},[',' M{1} ',']);
    idxM=find(idxM==idx);
    
    M{3} = DB_epm.data(:,idxM);
    
else

    for i=1:length(DB_epm.textdata(1,:))
      idx=strcmp(DB_epm.textdata{1,i},M{1});
      if idx==1
        M{3}=DB_epm.data(:,i-1);
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB_TB_V.textdata(1,:))==1
  
    DB_TB_V.textdata{1,1}=[DB_TB_V.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB_TB_V.textdata{1,1},',');
    idx=strfind(DB_TB_V.textdata{1,1},[',' M{1} ',']);
    idxM=find(idxM==idx);

    M{4} = DB_TB_V.data(:,idxM);

else

    for i=1:length(DB_TB_V.textdata(1,:))
      idx=strcmp(DB_TB_V.textdata{1,i},M{1});
      if idx==1
        M{4}=DB_TB_V.data(:,i-1);
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(DB_TB_K.textdata(1,:))==1

    DB_TB_K.textdata{1,1}=[DB_TB_K.textdata{1,1} ',']; % patch, add a comma "," at the end
    idxM=strfind(DB_TB_K.textdata{1,1},',');
    idx=strfind(DB_TB_K.textdata{1,1},[',' M{1} ',']);
    idxM=find(idxM==idx);
    
    M{5} = DB_TB_K.data(:,idxM);
    
else

    for i=1:length(DB_TB_K.textdata(1,:))
      idx=strcmp(DB_TB_K.textdata{1,i},M{1});
      if idx==1
        M{5}=DB_TB_K.data(:,i-1);
        % break % removing the break makes it slower but more compatible between Matlab and Octave
      end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%