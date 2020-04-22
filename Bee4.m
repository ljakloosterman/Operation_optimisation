clear all; close all; clc

%% Generate problem data
tic
rng(3,'twister')                            % random seed
D = 60;                                      % demand
range = 0.1;                                % range of facilities
ncities = 100;                               % # cities
nfacilities = 50;                            % # facilities

citiesx = zeros(ncities,1);                 % initialize x-cities
citiesy = citiesx;                          % initialize y-cities
facilitiesx = zeros(nfacilities,1);         % initialize x-facilities
facilitiesy = facilitiesx;                  % initialize y-facilities
rfacilities = range*ones(nfacilities,1);    % facilities radius

n=1;
while (n <= ncities)                
    citiesx(n) = rand*1.5;         
    citiesy(n) = rand;         
    n = n+1;     
end

depotx=rand*1.5;
depoty=rand;

n=1;
while (n <= nfacilities)              
    facilitiesx(n) = rand*1.5;         
    facilitiesy(n) = rand;         
    n = n+1;     
end

figure(1);                                    
plot(citiesx,citiesy,'*b')
hold on  
plot(facilitiesx,facilitiesy,'g.','MarkerSize',36)
plot(depotx,depoty,'r.','MarkerSize',36)
viscircles([facilitiesx,facilitiesy],rfacilities)

%% Compute objective function

nstops=nfacilities+1;
stopsx=[depotx;facilitiesx];    % X locations of depot and facilities
stopsy=[depoty;facilitiesy];    % Y locations of depot and facilities

idxs = nchoosek(1:nstops,2);    % All possible edges between stops

dist = hypot(stopsx(idxs(:,1)) - stopsx(idxs(:,2)), ... % Distances of al possible edges
             stopsy(idxs(:,1)) - stopsy(idxs(:,2)));
nedges = length(dist);                                   % Number of possible edges
nvariables=nedges+nfacilities;                           % Number of variables

f=[dist;zeros(nfacilities,1)];  % Adding a zero for each facility to obtain ojective function

%% Demand constraint (inequality)

% Add zero row to make room for demand constraints
Ain=zeros(1,nedges+nfacilities);

% Compute distances from city i to al facilities
for i=1:ncities
    dtab=zeros(1,nfacilities);
    for j=1:nfacilities
        dtab(j)=sqrt((citiesx(i)-facilitiesx(j))^2+(citiesy(i)-facilitiesy(j))^2); 
    end
    % Check whether closest facilities lies within range
    % If true add -1 to value with index of the facility
    if min(dtab)<range
        Ain(nedges+find(min(dtab)==dtab))=Ain(nedges+find(min(dtab)==dtab))-1;
    end
end

Bin=-D;

%% Equality constraints
% If a facility is visited two edges should be connected to it
% Two edges should be connected to the depot

% Add zero row to make room for demand constraints
Aeq=zeros(nstops,nvariables);

% Adding 1 for each edge which is connected to facility i
start=1;
len=nfacilities-1;
for i=1:nfacilities
    Aeq(i,start:start+len)=1;
    Aeq(i+1:nstops,start:start+len)=eye(len+1);
    start=start+len+1;
    len=len-1;
end

% Adding -2 to facility i on each row
n=nfacilities;
for i=1:nfacilities
    Aeq(i+1,nedges+i)=-2;
end

Beq=2;                          % Edges connected to depot is 2
Beq=[Beq;zeros(nfacilities,1)]; % If a facility is visited constraint is equal to 0     

%% Solving MILP

intcon = 1:nvariables;
lb = zeros(nvariables,1);   % Lower bound
ub = ones(nvariables,1);    % Upper bound

[xmin, fmin,exitflag]=intlinprog(f,intcon,Ain,Bin,Aeq,Beq,lb,ub);

%% Plot optimal path with subtours

lh = zeros(nstops,1);
lh = updateSalesmanPlot(lh,xmin(1:nedges),idxs,stopsx,stopsy);

%% Find subtours

% Detect subtours and store each tour in one array cell
tours = detectSubtours(xmin(1:nedges),idxs);
nsubtours = length(tours);                       % # if subtours
fprintf('# of subtours: %d\n',nsubtours);

%% Add subtour constraint and solve again
% A subtour will always have the same number of edges as facilities.
% Detecting a subtour and adding a constraint on all possible edges between
% these stops will eliminate the subtour. Between the subtour stops there
% should exist less edges than the number of stops.

% repeat until there is just one subtour
while nsubtours > 1 
    % Add zero rows to make room for subtour constraints
    Ain = [Ain;zeros(nsubtours,nvariables)];
    Bin = [Bin;zeros(nsubtours,1)];
    
    % Iterating through each subtour
    for i = 1:nsubtours
        newrow = size(Ain,1)+1;                     % Counter for indexing
        subtour = tours{i};                         % Extract the current subtour
        var = nchoosek(1:length(subtour),2);        % Find all possible connections between subtour facilities
        for j = 1:length(var)
            edge = (sum(idxs==subtour(var(j,1)),2)) & ... % Locating the subtour edges in list of all edges
                   (sum(idxs==subtour(var(j,2)),2));
            Ain(newrow,edge) = 1;                         % Adding 1 to constraint at edge locations
        end
        Bin(newrow) = length(subtour)-1;                  % One less trip than subtour stops
    end

    % Try to optimize again
    [xmin, fmin,exitflag]=intlinprog(f,intcon,Ain,Bin,Aeq,Beq,lb,ub);
    
    % Visualize result
    lh = updateSalesmanPlot(lh,xmin(1:nedges),idxs,stopsx,stopsy);
    
    % How many subtours this time?
    tours = detectSubtours(xmin(1:nedges),idxs);
    nsubtours = length(tours); % number of subtours
    fprintf('# of subtours: %d\n',nsubtours);
end
toc