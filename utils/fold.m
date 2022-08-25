function fold
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% This Matlab/Octave script is designed to fold a k-path into the BZ of 
% a supercell and produce KPOINTS file for VASP and case.klist_band file 
% for WIEN2k. The wavefunctions generated with VASP can be unfolded back 
% to a desired k-path set below.
%
% (c) Oleg Rubel, modified Jul 30, 2022
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clear all;

%% User input

kpath = [0 0 0
         0.5 0 0.25
         0.5 0.5 0.5
         0 0 0]; % desired k-path after unfolding
npath = [40 28 28]; % # either the total number of points along each segment
Dp2s = [1 -1 -2
        1 1 -2
        0 0 4]; % matrix used to transform a primitive cell to a supercell
toldk = 1e-6; % tolerance for round off errors in k values


%% Check input

if (size(kpath,1)-1 ~= length(npath))
    msg = MException('f2b:inpcheck',...
        'dimensions "kpath" and "npath" do not agree: %d vs %d',...
        size(kpath,1)-1, length(npath));
    throw(msg);
end
if det(Dp2s) <= 0 % Positive Definite Matrix
    msg = MException('f2b:inpcheck',...
        join(['Transformation matrix Dp2s is not a ',...
        'positive definite matrix: det(Dp2s) = %d']),...
        det(Dp2s));
    throw(msg);
end
if abs(det(Dp2s) - int16(det(Dp2s))) > toldk % Positive Definite Matrix
    msg = MException('f2b:inpcheck',...
        join(['Transformation matrix Dp2s does ',...
        'not give an integer volume scale: det(Dp2s) = %d']),...
        det(Dp2s));
    throw(msg);
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
% eliminate duplicates
if exist('OCTAVE_VERSION', 'builtin') ~= 0 % Octave
    kpr=unique(kpr,'rows'); % eliminate and sort
else % MATLAB
    kpr=unique(kpr,'rows','stable'); % no sorting
end
npt=size(kpr,1); % recalculate number of k-points

ksc = zeros(size(kpr)); % allocate k-list for supercell
for i=1:npt
    ksc(i,:) = kpr(i,:)*Dp2s;
end

% bring k-points coordinates to the range [0, 1)
for i=1:npt
    for j=1:3
        ksc(i,j) = mod( ksc(i,j), 1);
        if 1-ksc(i,j) < toldk
            ksc(i,j) = 0;
        end
    end
end

% clean up
ksc=round(ksc*1e14)/1e14; % round to 14 decimal points
% eliminate duplicates
if exist('OCTAVE_VERSION', 'builtin') ~= 0 % Octave
    ksc=unique(ksc,'rows'); % eliminate and sort
else % MATLAB
    ksc=unique(ksc,'rows','stable'); % no sorting
end
npt=size(ksc,1); % recalculate number of k-points
msg = [num2str(npt),' unique k points were generated.'];
disp(msg);


%% Write VASP KPOINTS file

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
disp('Generated k points are stored in KPOINTS file (for VASP)');


%% Write WIEN2k case.klist_band file

fileID = fopen('case.klist_band','w');
for i=1:npt
    % convert real coordinate of k points the ratio of integers
    ndrat = real2rat(ksc(i,:), 9); % 9 is related to the format output i9
    format = '%10i %9i %9i %9i %9i %4.2f\n';
    fprintf(fileID, format, i, ndrat, 1.0);
end
fprintf(fileID, '%s', 'END'); % append the k-list file
fclose(fileID);
disp('Generated k points are stored in case.klist_band file (for WIEN2k)');

% Nested functions
% -------------------------------------------------------------------------
    function out = real2rat(V, outsizemax)
    % transform vector V(1:3) into the ratio of integers nout(1:3)/dout
        [n,d] = rat(V); % bring rational numbers to the form n/d
        dout = lcm(d(1),d(2)); % determine a common multiplier
        dout = lcm(d(3),dout);
        nout = n*dout./d; % scale by the common multiplier
        out = [nout, dout]; % append common divisor
        outsize = max(ceil(log10(abs(out))),1); % size of the output integer
        % check if the integers are compatible with the output format
        if any(outsize > outsizemax) 
            disp('out ='); disp(outsize);
            disp('outsize ='); disp(outsize);
            disp('outsizemax ='); disp(outsizemax);
            msg = 'Integer size exceeds the output format.';
            error(msg);
        end
        % check accuracy of the transform
        errorOK = 1e-6;
        errork = V - nout./dout;
        if any(abs(errork) > errorOK)
            disp('WARNING:')
            disp('  k ='); disp(V);
            disp('  k1 k2 k3 kdivisor ='); disp(out);
            disp('  k1/kdivisor k2/kdivisor k3/kdivisor ='); disp(nout./dout);
            disp('  max abs k error ='); disp(abs(errork));
            disp('  tolerable error ='); disp(errorOK);
        end
    end % function real2rat
end % function fold