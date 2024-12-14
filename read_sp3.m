function sp3=read_sp3(fname)
% function data=readsp3_pa(fname)
% P. Axelrad
% Y. Wang edited for multi-constellation (2023-09)
%
% Function to read sp3 files
% Input is the filename
%
% Output is array of the data 
% Column 
%   1 - Week number
%   2 - TOW (s)
%   3 - PRN
%   4-6 - XYZ Position (km)
%   7 - Clock bias (microsec)
%   8 - constellation (1-GPS, 2-GLO, 3-GAL, 4-BDS, 5-QZSS)
% Current version does not output any of the accuracy values.
% Functions called - cal2gps (by P. Axelrad)
%

sp3=[];
fin = fopen(fname,'r');

line = fgets(fin);
while ( line ~= -1 )
    msg=line(1);
    switch ( msg )
        case '*'
            temp=sscanf(line(4:end),'%f')';
            [wn,t]=cal2gps(datetime(temp));
        case 'P'
            constellation = find(line(2)==['G','R','E','C','J']);
            prn_num = str2num(line(3:4));
            temp=sscanf(line(6:end),'%f')';
           	sp3=[sp3;[wn t prn_num temp(1:4) constellation]];
    end
    line = fgets(fin);
end
return
