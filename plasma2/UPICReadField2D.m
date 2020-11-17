%clear all

%Basic plotting diagnostic for 2D UPIC
%parameters from 2D UPIC input file
indx = 9;
indy = 3;
tend = 40.08;
dt = 0.08;
nt = 10;
np = 1;
dx = 1;
dy = 1;

%options
read = 1;
read_start = 1;
read_stop = 501;


if read == 1
steps = round((read_start-read_stop)/nt) + 1;

  disp('read in files from 2-D UPIC')

  %read in field data from bbeps2
  disp(' -field data')
  %ex_txy = zeros(2^indx,2^indy,steps);
  %ey_txy = zeros(2^indx,2^indy,steps);
%  ez_txy = zeros(2^indx,2^indy,steps);
%  bx_txy = zeros(2^indx,2^indy,steps);
%  by_txy = zeros(2^indx,2^indy,steps);
  %bz_txy = zeros(2^indx,2^indy,steps);
%  energy = zeros(steps);

  for i = read_start:nt:read_stop

    step = i-1

%    ex_txy(:,:,step+1) = emma_2d('MS','e1',1,step,indx,indy,2);
%    ey_txy(:,:,step+1) = emma_2d('MS','e2',1,step,indx,indy,2);
    ez_txy(:,:,step+1) = emma_2d('MS','e3',1,step,indx,indy,2);

    %bx_txy(:,:,step+1) = emma_2d('MS','b1',1,step,indx,indy,2);
%    by_txy(:,:,step+1) = emma_2d('MS','b2',1,step,indx,indy,2);
%    bz_txy(:,:,step+1) = emma_2d('MS','b3',1,step,indx,indy,2);

  end

end