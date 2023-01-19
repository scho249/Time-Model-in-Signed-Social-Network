clear; close all; clc;
%% triad


X0 = [0 1 1 0;
      1 0 -1 0;
      1 -1 0 0;
      0 0 0 0];
 X0 = triu(X0); sz = size(X0(:,1));
 pos_edge = 0;
 neg_edge = 0;
 for i = 1:sz
     if pos_edge + neg_edge < 3
        for j = 1:sz  
            if X0(i,j) > 0
                 pos_edge = pos_edge + X0(i,j);
             elseif X0(i,j) < 0
                neg_edge = neg_edge + X0(i,j);
            end
        end
     end
 end
 i3 = (pos_edge*neg_edge)/(pos_edge+abs(neg_edge));
 