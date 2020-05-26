function ubs_dots
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot undolded band structure with matlab using boubles
% 
% (c) Oleg Rubel, modified May 26, 2020
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Init. parameters
KPATH = [0 0 0;...
         1/2 1/2 1/2; ...
         1/2 1/2 0]; % k-point path
FOLDS = [4 4 4]; % multiplicity in the corresponding directions used when constructing the super-cell
KLABEL = {'G'; 'R'; 'S'};
finpt = 'WAVECAR_spin1.f2b'; % input file name
Ef = 1.325692; % Fermi energy (Ry)
ERANGE = [Ef-1 Ef+2]; % energy range for plot (eV)
pwr = 1/1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.
         % 0 - folded bands (needs wth = 0)
msz = 10; % marker size for plot
lwdth = 0.5; % plot line width
fontSize = 9; % points
PLTSZ = [1 1 600/1.5 300/1.5]; % plot size
wth = 0.05; % threshold weight
clrmp = jet;    % flipud(gray)
                % flipud(pink)
                % flipud(hot)
                % flipud(autumn)
                % cool
                % flipud(bone)
                % flipud(jet)
                % jet
G = [0.039913469  0.000000000  0.000000000
0.000000000  0.039913457  0.000000000
0.000000000  0.000000000  0.039450845];    % Reciprocal latt. vect. from *.outputkgen


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
XTICKS = [0];
for ikp = 1 : size(KPATH,1)-1
    B = KPATH(ikp,:) - KPATH(ikp+1,:);
    dk = sqrt(dot(B,B));
    XTICKS = [XTICKS; XTICKS(ikp)+dk];
    for j = 1 : length(EIG)
        if EIG(j) > ERANGE(1) && EIG(j) < ERANGE(2) && W(j) >= wth
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
            if dist < eps % k-point is on the path
                A = KPATH(ikp,:) - KEIG(j,:);
                x = dot(A,B)/dk;
                if x >= 0  &&  x <= dk+eps % k-point is within the path range
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


%% Plot results
hFig = figure(1);

% Fig 1(a)
subplot(1,2,1);
set(gca,'FontSize',fontSize);
set(hFig, 'Position', PLTSZ, 'PaperPositionMode','auto')
map = colormap(clrmp);
WGHTRS = rescale(WGHT,pwr);
scatter(L,(ENE-Ef), WGHTRS*msz, WGHTRS,'LineWidth',lwdth);
hold on;
axis([XTICKS(1) XTICKS(end) ERANGE(1)-Ef ERANGE(2)-Ef])
yticks = get(gca,'ytick');
set(gca,'YTick',yticks);
for i = 1 : length(yticks)
    newYTick{i} = sprintf('%1.1f',yticks(i));
end
set(gca,'YTickLabel',newYTick);
hline = plot([0 XTICKS(end)],[0 0]); % Fermi level
set(hline,'Color','k','LineStyle','--');
set(gca,'XTick',XTICKS);
set(gca,'XTickLabel',KLABEL);
set(gca,'XGrid','on', 'GridLineStyle','-');
caxis([0 1]); % normalize intensities to 1
xlabel('Wave vector')
ylabel('Energy (eV)')
box on
hold off

% Fig 1(b)
subplot(1,2,2);
set(gca,'FontSize',fontSize);
DAT = linspace(0,1,10);
DATX = ones(size(DAT));
DATRS = rescale(DAT,pwr);
scatter(DATX,DAT, DATRS*msz, DATRS,'LineWidth',lwdth);
caxis([0 1])
ylabel('Spectral weight')

% SAVE plot as *.eps
print( [finpt '.eps'], '-depsc')

% -------------------------------------------------------------------------
function W = coordTransform(V,G)
% transform vector V(:,3) in G(3,3) coord. system -> W(:,3) in Cartesian coordinates
% G vector elements are in columns!
W = zeros(size(V));
for i = 1:size(V,1)
    W(i,:) = G(1,:)*V(i,1) + G(2,:)*V(i,2) + G(3,:)*V(i,3);
end;
% -------------------------------------------------------------------------
function WRESCL = rescale(W,pwr)
% rescale weights using a power functio W^pwr
WRESCL=W.^(pwr); % rescale if needed to enhance
WRESCL = WRESCL + eps; % need eps to make plot "heapy"
% -------------------------------------------------------------------------
function [KEIG, EIG, W] = readinput(filename)
% read input data
DATA = importdata(filename);
KEIG = DATA(:,1:3);
EIG = DATA(:,4);
W = DATA(:,5);
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
% -------------------------------------------------------------------------
