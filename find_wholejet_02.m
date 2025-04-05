function [wholejet_value,wholejet_x,wholejet_y,jetaxis_value,jetaxis_x,jetaxis_y] = ...
    find_wholejet_01(inputdata,inputx,inputy,jet_distinguish,jet_length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the jets axes in two dimensional matrix and separate them
% -------------------------------------------------------------------------
% input:
%       inputdata: 2D matrix
%       inputx:    longitude of the matrix
%       inputy:    latitude of the matrix
%       jet_distinguish: separate two jets whose distance is large than jet_distinguish
%       jet_length: save the jets whose length is less than jet_length
%
% output:
%       wholejet_value (jetnumber*length*width):all jets value include
%       length along the axis and width.
%       wholejet_x (jetnumber*length*width): x coordinate
%       wholejet_y (jetnumber*length*width): y coordinate
%       jetaxis_value (jetnumber*length): the matrix contain all the jets found in 'inputdata'
%       jetaxis_x (jetnumber*length):     x coordinate of the 'jetaxis_value'
%       jetaxis_y (jetnumber*length):     y coordinate of the 'jetaxis_value'
% -------------------------------------------------------------------------
% by zhangyu 20200810
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====================== find extreme points ============================ %
inputdata(isnan(inputdata)==1) = 0; 
extreme_point = zeros(size(inputdata));
extreme_y = zeros(size(inputdata));
extreme_x = zeros(size(inputdata));
for i = 1:length(inputx)
    tmp = [inputdata(:,i);0];
    k = 0;
    for j = 1:length(tmp)
        if tmp(j) ~= 0 && tmp(j+1) ~= 0
            tmp_new(j,1) = tmp(j);
            tmp(j) = 0;

        elseif tmp(j) ~= 0 && tmp(j+1) == 0
            tmp_new(j,1) = tmp(j);
            tmp(j) = 0;
            
            k = k + 1;      
            % find the extreme points
            tmp_new(tmp_new==0) = [];  % non-zero values
            if tmp_new > 0 
                ex_point = max(tmp_new);
            elseif tmp_new < 0
                ex_point = min(tmp_new);
            end
            extreme_point(k,i) = ex_point;
            ex_y = inputy( find( abs(inputdata(:,i) - ex_point) < 1.e-5 ));extreme_y(k,i) = ex_y(1);
            extreme_x(k,i) = inputx(i);
            
            clear tmp_new
        end
    end        
end

% ======================= find jets axis ================================ %
jetaxis_x = zeros(size(inputdata));
jetaxis_y = zeros(size(inputdata));
jetaxis_value = zeros(size(inputdata));
jet_num = 0;
extreme_point_new = extreme_point;
%size(extreme_point_new)
%size(inputy)

for i = 1:length(inputx)-1
    for j = 1:length(inputy)
        
        % find the start point of the jet ----------------------- %
        if extreme_point_new(j,i) ~= 0 
            jet_num = jet_num + 1;   % jet number
            
            % detect the jet axes from the start point
            m = j; n = i;jet_len = 1;
            for t = 1:length(inputx) % iterative in dimention t whose length is same as inputx  
                for k = 1:length(inputy)
                    if extreme_point_new(m,n) ~= 0 && extreme_point_new(k,n+1) ~= 0 ...
                            && abs(extreme_y(m,n) - extreme_y(k,n+1)) <= jet_distinguish
                        jetaxis_value(jet_num,jet_len) = extreme_point_new(m,n);
                        jetaxis_y(jet_num,jet_len) = extreme_y(m,n);
                        jetaxis_x(jet_num,jet_len) = extreme_x(m,n);
                        extreme_point_new(m,n) = 0;
                        
                        if n < length(inputx)-1
                            m = k; 
                            n = n+1;
                            jet_len = jet_len + 1;   % jet length
                        end
                        
                    end
                end
            end
            % -------------------------------------- %
            
        end
        % ------------------------------------------------------ %
        
    end
end

% choose the jets by control jet length
index_y = zeros(size(inputy));
for le = 1:size(jetaxis_value,1)
    if mean(jetaxis_value(le,:)) == 0  % jet length larger than 400 km
        index_y(le) = 1;
    else
        tmp = jetaxis_value(le,:);tmp(tmp==0) = [];
        if length(tmp) < jet_length    % jet length larger than 400 km
            index_y(le) = 1;
        end
    end
    clear tmp
end

jetaxis_value(find(abs(index_y-1) < 1.e-5),:) = [];
jetaxis_x(find(abs(index_y-1) < 1.e-5),:) = [];
jetaxis_y(find(abs(index_y-1) < 1.e-5),:) = [];

% ====================== find whole jets ================================ %
wholejet_value = zeros(size(jetaxis_value,1),size(jetaxis_value,2),length(inputy));
wholejet_x = zeros(size(jetaxis_value,1),size(jetaxis_value,2),length(inputy));
wholejet_y = zeros(size(jetaxis_value,1),size(jetaxis_value,2),length(inputy));

inputdata1 = [inputdata;zeros(1,size(inputdata,2))];
for i = 1:size(jetaxis_value,1)
    for j = 1:size(jetaxis_value,2)
%         disp(['i=',num2str(i),'j=',num2str(j)])
        
       if jetaxis_value(i,j) ~= 0
           m = find(abs(inputy - jetaxis_y(i,j))<1.e-5);  % jetaxis point location
           n = find(abs(inputx - jetaxis_x(i,j))<1.e-5);
           
           jetwidth_value = zeros(size(inputy));
           jetwidth_y = zeros(size(inputy));
           
           % search downward
           for t = 1:length(inputy)
               if m > 1
                   if inputdata1(m,n) ~= 0 && inputdata1(m-1,n) ~= 0
                       jetwidth_value(m) = inputdata(m,n);
                       jetwidth_y(m) = inputy(m);
                   elseif inputdata1(m,n) ~= 0 && inputdata1(m-1,n) == 0
                       jetwidth_value(m) = inputdata(m,n);
                       jetwidth_y(m) = inputy(m);
                       break
                   end
                   m = m - 1;
               end
           end
           
           % search upward
           for t = 1:length(inputy)
               if m < length(inputy)-1
                   if inputdata1(m,n) ~= 0 && inputdata1(m+1,n) ~= 0
                       jetwidth_value(m) = inputdata(m,n);
                       jetwidth_y(m) = inputy(m);
                   elseif inputdata1(m,n) ~= 0 && inputdata1(m+1,n) == 0
                       jetwidth_value(m) = inputdata(m,n);
                       jetwidth_y(m) = inputy(m);
                       break
                   end
                   m = m + 1;
               end
           end
           
           jetwidth_x = inputx(n)*ones(size(jetwidth_y));
           
           wholejet_value(i,j,:) = jetwidth_value;
           wholejet_x(i,j,:) = jetwidth_x;
           wholejet_y(i,j,:) = jetwidth_y;
           
           clear jetwidth_x jetwidth_y jetwidth_value

       end 
       
      
    end
end
end
