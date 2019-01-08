% This program computes the Reciprocal lattice vectors Gvec that is needed in the 
% epm model. According to Dragica Vasileska (Arizona Satae University), the reciprocal
% lattice vector up to and includingthe 10th-nearest neighbours are sufficient and correspond 
% to 137 plane waves
% https://nanohub.org/resources/1524/download/empiricalpseudopotentialmethod_word.pdf


close all
clear all
clc

NN_neighbours=10;
ig2_mag=[0, 3, 4, 8, 11, 12, 16, 19, 20, 24];
ig2_mag=ig2_mag(1:NN_neighbours);

icount=0;

for i=1:NN_neighbours
ig2=ig2_mag(i);
ia=floor(sqrt(double(ig2)));

ic=0;
for i=-ia:ia
  for j=-ia:ia
    for k=-ia:ia
         
      idx = i^2 + j^2 + k^2;
       
      if (idx==ig2) 
        ic=ic+1;
        gx(ic+icount)=i;
        gy(ic+icount)=j;
        gz(ic+icount)=k;
      end
      
      
    end
  end
end

icount=icount+ic;

end

Gvec=[gx' gy' gz'];

size(Gvec)
save('Gvec.txt','Gvec')