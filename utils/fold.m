%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% This Matlab script is designed to fold a k-path into the BZ of 
% a supercell and produce KPOINTS file for VASP. The wavefunctions 
% generated with VASP can be unfolded back to a desured k-path set 
% below.
%
% (c) Oleg Rubel, modified Oct 7, 2017
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


clear all;

%% User input
kpath = [1/2 0 0; ...
         0 0 0;...
         1/2 1/2 0];
npath = [10 15]; % # of points along each segment
folds = [1 2 3]; % multiplicity used to create a supercell

%% Check
if (size(kpath,1)-1 ~= length(npath))
    error('dimensions "kpath" and "npt" do not agree')
end

%% Compute
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

figure(1);
plot(kpr(:,1),kpr(:,2),'o')
grid on
pbaspect([1 1 1])

ksc = zeros(size(kpr)); % allocate k-list in supercell
for i=1:npt
    ksc(i,:) = kpr(i,:).*folds;
end

figure(2);
plot(ksc(:,1),ksc(:,2),'o')
grid on
pbaspect([1 1 1])

% bring the result between [-0.5, 0.5]
for i=1:npt
    for j=1:3
        a = ksc(i,j);
        b = a+1/2;
        c = floor(b);
        ksc(i,j) = ksc(i,j) - c;
    end
end

figure(3);
plot(ksc(:,1),ksc(:,2),'o')
grid on
pbaspect([1 1 1])

%% Write KPOINTS file

fileID = fopen('KPOINTS','w');
fprintf(fileID,'%34s\n','k-mesh for unfolder band structure');
fprintf(fileID,'%8i\n',npt);
fprintf(fileID,'%18s\n','Reciprocal lattice');
for i=1:npt
    if (i ~= npt)
        fprintf(fileID,'%13.10f %13.10f %13.10f  1\n',ksc(i,:));
    else
        fprintf(fileID,'%13.10f %13.10f %13.10f  1',ksc(i,:));
    end
end
fclose(fileID);
