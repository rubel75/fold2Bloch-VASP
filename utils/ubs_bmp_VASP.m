function ubs_bmp_VASP
% Prepare a binary file to plot the undolded band structure as a bitmap figure
%
% Execution:
%    octave ubs_bmp_VASP
%
% Generated file:
%    WAVECAR_spin1.f2b.bin
%
% Band structure plotting:
%    gnuplot f2b-band-structure.plt
%
% (c) Oleg Rubel modified 26 May 2020

%% Init. parameters
KPATH = [0 0 1/2;...
         0 0 0; ...
         1/2 1/2 0]; % k-point path
FOLDS = [1 1 6]; % multiplicity in the corresponding directions used when constructing the super-cell
KLABEL = {'Z'; 'G'; 'L'};
finpt = 'WAVECAR_spin1.f2b'; % input file name
Ef = 3.237409; % Fermi energy (Ry)
ERANGE = [Ef-2.0 Ef+3.0]; % energy range for plot (eV)
pwr = 1/1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.
         % 0 - folded bands (needs wth = 0)
nK = 280; % pixels along k
nE =  400; % pixels along Energy axis
sK = 0.0025; % smearing factor in k-space
sE = 0.025; % smearing factor in energy
G = [0.259148502 -0.149619435 -0.104836011
0.000000000  0.299238904 -0.104836011
0.000000000  0.000000000  0.052418005];    % Reciprocal latt. vect. from OUTCAR
roundOffErrK = 0.000001; % this is the round off error 1/3 = 0.333333 + err

%% INITIALIZATION
[KEIG, EIG, W] = readinput(finpt); % read input data from file
% EIG - energy eigenvalues
% KEIG - k-list for eigenvalues
% W - list of characters

%% MAIN
L = [];
ENE = [];
WGHT = [];
%G = G'; % transpose G matrix (need for Wien2k)
for i=1 : 3
    G(i,:)=G(i,:)*FOLDS(i); % rescale reciprocal lattice vectors 
end                         % from supercell to primitive cell
dl = 0; % cumulative length of the path
KPATH = coordTransform(KPATH,G);
KEIG = coordTransform(KEIG,G);
epsk = [roundOffErrK roundOffErrK roundOffErrK]; % k rounding error
epsk = coordTransform(epsk,G); % transform to Cart. coords
epsk = sqrt(dot(epsk,epsk)); % get magnitude of the vector
XTICKS = [0];
for ikp = 1 : size(KPATH,1)-1
    B = KPATH(ikp,:) - KPATH(ikp+1,:);
    dk = sqrt(dot(B,B));
    XTICKS = [XTICKS; XTICKS(ikp)+dk];
    for j = 1 : length(EIG)
        if EIG(j) > ERANGE(1) && EIG(j) < ERANGE(2) && W(j) > 0
            dist = Inf; % initialize distance to a path
            for ikx = -1:1 % include periodic images of the BZ
                for iky = -1:1
                    for ikz = -1:1
                        KPERIOD = [ikx iky ikz]; % periodic shift
                        % transform to Cartezian coords
                        KPERIOD = coordTransform(KPERIOD,G);
                        % evaluate distance to the path
                        dist2 = dp2l( KEIG(j,:) + KPERIOD , ...
                            KPATH(ikp,:) , KPATH(ikp+1,:) );
                        % select smallest distance
                        dist = min(dist,dist2);
                    end
                end
            end
            if dist < epsk % k-point is on the path
                A = KPATH(ikp,:) - KEIG(j,:);
                x = dot(A,B)/dk;
                if x >= 0  &&  x <= dk+epsk % k-point is within the path range
                    L = [L; x+dl]; % append k-point coordinate along the path
                    ENE = [ENE; EIG(j)]; % append energy list
                    WGHT = [WGHT; W(j)];
                end
            end
        end
    end
    dl = dl + dk;
end
if isempty(L)
    msg = ['No eigenvalues are selected for the plot. ', ...
        'The likely reason is that the energy range is ', ...
        'too restrictive (check ERANGE), or no k-points are located ', ...
        'on the path selected (check KPATH)'];
    error(msg);
end


WGHTRS = rescale(WGHT,pwr);
Li = linspace(XTICKS(1),XTICKS(end),240);
Ei = linspace(ERANGE(1)-Ef,ERANGE(2)-Ef,200);
nLi = length(Li);
nEi = length(Ei);
Wi = zeros(nLi,nEi);
ENE = ENE-Ef;
[XTICKS(1) XTICKS(end)]
for l = 1:length(L) % loop over all unfolded points
    disp(100*l/length(L));
    for i = 1:nLi
    if abs(Li(i)-L(l)) > 3*sK
        continue; % skip if > 3*sigma
    endif
        for j = 1:nEi
            if abs(Ei(j)-ENE(l)) > 3*sE
                continue; % skip if > 3*sigma
            endif
            wK = gauss(Li(i),L(l),sK);
            wE = gauss(Ei(j),ENE(l),sE);
            Wi(i,j) = Wi(i,j) + wK*wE*WGHTRS(l);
        end
    end
end
Wi = Wi/(gauss(0,0,sK)*gauss(0,0,sE)); % normalize weights
save_binary_matrix([finpt,'.bin'],Li,Ei,Wi')
XTICKS
KLABEL

%set(gca,'XTick',XTICKS);
%set(gca,'XTickLabel',KLABEL);

% -------------------------------------------------------------------------
function RES = gauss(X,x0,s) % bound. conditions
RES = 1/(s*sqrt(2*pi)) * exp(-((X-x0).^2)/(2*s^2));
endfunction
% -------------------------------------------------------------------------
function W = coordTransform(V,G)
% transform vector V(:,3) in G(3,3) coord. system -> W(:,3) in Cartesian coordinates
% G vector elements are in columns!
W = zeros(size(V));
for i = 1:size(V,1)
    W(i,:) = G(1,:)*V(i,1) + G(2,:)*V(i,2) + G(3,:)*V(i,3);
end;
end
% -------------------------------------------------------------------------
function WRESCL = rescale(W,pwr)
% rescale weights using a power functio W^pwr
WRESCL=W.^(pwr); % rescale if needed to enhance
WRESCL = WRESCL + eps; % need eps to make plot "heapy"
end
% -------------------------------------------------------------------------
function [KEIG, EIG, W] = readinput(filename)
% read input data
DATA = importdata(filename);
KEIG = DATA(:,1:3);
EIG = DATA(:,4);
W = DATA(:,5);
end
% -------------------------------------------------------------------------
function RES = dp2l(X0,X1,X2) % distance from point {X0} to line {X1}-{X2}
% see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
denom = X2 - X1;
denomabs = sqrt(dot(denom,denom));
if denomabs < eps
    display(X1); display(X2);
    error('X1 = X2');
end;
numer = cross( X0-X1 , X0-X2 );
numerabs = sqrt(dot(numer,numer));
RES = numerabs/denomabs;
end
% -------------------------------------------------------------------------
function save_binary_matrix(file,x,y,M)
% SAVE_BINARY_MATRIX save x,y,M in the matrix binary format of Gnuplot
%   Usage: save_binary_matrix(file,x,y,M)
%
%   Input parameters:
%       file    - filename of the data file
%       x       - x axis values
%       y       - y axis values
%       M       - matrix data size(M)=y,x
%
%   SAVE_BINARY_MATRIX(file,x,y,M) saves the values of x,y and M in a binary
%   matrix format useable by Gnuplot, see
%   http://www.gnuplot.info/docs_4.4/gnuplot.html#x1-33700077.1.1 for details.
%

% AUTHOR: Hagen Wierstorf


% Checking of input  parameters 
error(nargchk(4,4,nargin));


% Computation 

% Check if the data has the right format
[ly,lx] = size(M);
if lx~=length(x) || ly~=length(y)
    error('%s: size(M) has to be y,x!',upper(mfilename));
end

% Create matrix to store in the file
MS = zeros(length(x)+1,length(y)+1);
MS(1,1) = length(x);
MS(1,2:end) = y;
MS(2:end,1) = x;
MS(2:end,2:end) = M';

% Write data into the file
fid = fopen(file,'w');
fwrite(fid,MS,'float');
fclose(fid);
end

end % function ubs_bmp_VASP
