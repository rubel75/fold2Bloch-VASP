%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This Matlab script is designed to fold a k-path into the BZ of 
% a supercell and produce KPOINTS file for VASP. The wavefunctions 
% generated with VASP can be unfolded back to a desired k-path set 
% below.
%
% (c) Oleg Rubel, modified Dec 06, 2018
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all;

%% User input

kpath = [0 0 0; ...
         1/2 1/2 1/2;...
         1/2 1/2 0]; % desired k-path after unfolding
npath = [17 10]; % # of points along each segment
folds = [2 2 2]; % multiplicity used to create a supercell


%% Check input

if (size(kpath,1)-1 ~= length(npath))
    error('dimensions "kpath" and "npath" do not agree')
end


%% Compute folded k-points

npt = sum(npath) - length(npath) + 1;
kpr = zeros(npt,3);
for i=1:size(kpath,1)-1
    kpr1 = kpath(i,:);
    kpr2 = kpath(i+1,:);
    dk = (kpr2 - kpr1)./(npath(i)-1);
    if (i > 1)
        index = sum(npath(1:i-1));
    else
        index = 1;
    end
    kpr(index,:) = kpr1;
    for j=index+1:index+npath(i)-1
        kpr(j,:) = kpr(j-1,:) + dk;
    end
end
kpr=round(kpr*1e14)/1e14; % round to 14 decimal points
kpr=unique(kpr,'rows','stable'); % eliminate duplicates
npt=size(kpr,1); % recalculate number of k-points

ksc = zeros(size(kpr)); % allocate k-list for supercell
for i=1:npt
    ksc(i,:) = kpr(i,:).*folds;
end

% bring k-points coordinates to the range [-0.5, 0.5]
for i=1:npt
    for j=1:3
        a = ksc(i,j);
        if a < 0
          b = a+1/2;
          c = floor(b);
        else
          b = a-1/2;
          c = ceil(b);
        end
        ksc(i,j) = a - c;
    end
end

% clean up
ksc=round(ksc*1e14)/1e14; % round to 14 decimal points
ksc=unique(ksc,'rows','stable'); % eliminate duplicates
npt=size(ksc,1); % recalculate number of k-points


%% Write KPOINTS file

fileID = fopen('KPOINTS','w');
fprintf(fileID,'%34s\n','k-mesh for unfolding the band structure');
fprintf(fileID,'%8i\n',npt);
fprintf(fileID,'%18s\n','Reciprocal lattice');
for i=1:npt
    if (i ~= npt) % print with return at the end "\n"
        fprintf(fileID,'%17.14f %17.14f %17.14f  0\n',ksc(i,:));
    else % no return for the last line
        fprintf(fileID,'%17.14f %17.14f %17.14f  0',ksc(i,:));
    end
end
fclose(fileID);
