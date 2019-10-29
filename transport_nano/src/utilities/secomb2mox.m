%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	MATLAB SCRIPT TO CONVERT THE SECOMB FORMAT OF THE NETWORK TO MOX FORMAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
T=rattm93b;
br=height(T);

fileID = fopen('rattm93b_coarse.pts','w');
fprintf(fileID,'%10s\r\n','BEGIN_LIST');
%writing in the same order of .pts branch 1 -> branch 23 , branch 23 ->
%branch 1, the rest equal
counter = 1;
for (b=[23, 2:22 , 1, 24:br])
    fprintf(fileID,'%9s\r\n','BEGIN_ARC');
    fprintf(fileID,'%3s\r\n','BC ');
    fprintf(fileID,'%3s\r\n','BC ');
    n = 11;
    x=(linspace(T{b,6}/50, T{b,9}/50, n));
    y=(linspace(T{b,7}/50, T{b,10}/50, n));
    z=(linspace(T{b,8}/50, T{b,11}/50, n));
    nb= counter*ones(1,n);
    fprintf(fileID,' %d \t  %4.8f \t %4.8f \t %4.8f \t start \r  ',[nb(1); x(1); y(1); z(1)]);
    fprintf(fileID,' %d \t  %4.8f \t %4.8f \t %4.8f \t end \r  ',[nb(end); x(end); y(end); z(end)]);
    fprintf(fileID,' %d \t  %4.8f \t %4.8f \t %4.8f \t point \r  ',[nb(2:end-1); x(2:end-1); y(2:end-1); z(2:end-1)]);
    fprintf(fileID,'%s\r\n','END_ARC');
    fprintf(fileID,'\r\n');
    fprintf(fileID,'\r\n');
    counter = counter +1;
end
fprintf(fileID,'%8s\r\n','END_LIST');
fclose(fileID);


fileID = fopen('rattm93a_norm50_radius.pts','w');
fprintf(fileID,'%10s\r\n','BEGIN_LIST');
r=T{:,4}/50./2.;
fprintf(fileID,' %4.8f \r  ',r);
fprintf(fileID,'%8s\r\n','END_LIST');

