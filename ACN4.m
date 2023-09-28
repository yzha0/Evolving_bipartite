%ACN4.m
%
%Idea is to start with a plant and animal mutualism and let each mutualism
%characteristic be on a 1D line.  Then there is some probability of
%speciation at each generation, and also some drift that occurs within each
%generation.  At each update connections are maintained, created and lost
%with some probability proportional to the trait distance between the
%mutualistic candidate partners. 
 
% edits made are that this does not keep track of more than 2 gens at a time
% and I am changing the matching function to a square instead of absolute value.
% I deleted cascade speciation
% I am now keeping track of events by adding a B matrix
% 7_22_11: changed selection differential
%with changes to fitness function
% 3_15_ changes visualization to make nodes appear in order of speciation
% Making corrections (changed alpha and beta) to reflect new KLS model,
% option 1, you look around and choose based on distribution of who is out
% there, not picking an interactor first.
% Asymmetric

%Final model

clear all;
clear figure;
clf;

extinctions = 0;
toobig = 0;

%parameters
T = 400; %stop time T
% think of the X species as fig
pX = 0.01; % probability of speciation, p_sx
% think of the Y species as fig wasp
pY = 0.01; %p_sy

maxNX =500;  % the maximum number of species for preallocation
maxNY =1500;

v = 1; % tolerance param for trait matching v2 in text, v^2
theta = .1; %assume thetaxi=thetayj=thetaxJ=thetayI
thetax =theta;%.1; %Vp phenotypic variance for x, theta^2_x
thetay =theta;%.1; %Vp phenotypic variance for y
xix = .1; %standard deviation for drift from one gen to the next fig 
xiy = .1; %standard deviation for drift from one gen to the next wasp
xis = .5; %standard deviation of trait value change for a speciation event
c = 0.99; %proportion of connection probability based on historical association,
         %so that if you were connected in the 
         %prior generation, then you will have a connection in this
         %generation
% alpha = v/(v+thetax); beta = v/(v+thetay); alpha is difined later 
% since it depends on which plant and which pollinator we are considering
         

%test null model where there is no trait variance either in drift or
%inheritance, all uncommented, total null, stays at 0, alpha=beta=0, then
%just random walk of unconnected species, which may eventually go extinct
%if c is not 1
%alpha = 0;
%beta = 0;
%xix=0;
%xiy=0;
%xis=0;
%c=1;

% x and y preallocated
x = NaN(2,maxNX); %preallocate matrix for speed
xall = NaN(T,maxNX);
y = NaN(2,maxNY);
yall = NaN(T,maxNY);

%Initial conditions for one species in each
x0 = 0; %inital trait value x of species X
x(1,1) = x0;
NsppX = 1; % total # of species of type X that have ever existed
iXextant = [1];  % the indicies of the extant species (out of the NsppX in x)
y0 = 0; %inital trait value y of species type Y
y(1,1)= y0;
NsppY = 1; % total # of species of type Y that have ever existed
iYextant = [1]; % the indicies of the extant species (out of the NsppyallY in y)
%B = zeros(maxNX,maxNY);
% We start with one species each that are connected to each other
A(1,1) = 1;
%B(1:1) = 1;
NextantsppX = length(iXextant);
NextantsppY = length(iYextant);

%Initial conditions for two species in each
%x(1,1:2) = [-1 1];
%NsppX = 2; % total # of species of type X that have ever existed
%iXextant = [1 2];  % the indicies of the extant species (out of the NsppX in x)
%y(1,1:2)= [-1 1];
%NsppY = 2; % total # of species of type Y that have ever existed
%iYextant = [1 2]; % the indicies of the extant species (out of the NsppyallY in y)
% We start out with two species connected to each other
%A(1:2,1:2) = [1 0;0 1];
%B(1:2,1:2)=[1 0;0 1];
%NoldsppXB = 1;
%NoldsppYB = 1;

%This recovery term would only kick in if I wanted to delay extinction
%after the loss of a connection
%recovery = 0;

% Keeping track of events
%cospeciation = 0;
%duplication = 0;
%host_switch = 0;
%sorting = 0;

%To keep track of phylogeny, we use an indexing pointer-like system  The
%first row holds the tree topology and the second row holds the metric of
%the tree (what generation speciation occured)
% For example, [1 3 2; 0 150 150] means
%       1     at start       
%      / \
%     /   \
%    1     2   at generation 100
%   / \     \
%  1   3     2 at generation 150

phylogenyX = NaN(2,NsppX);
%phylogenyX(:,1) = [1,0];
phylogenyX(:,1:2) = [1 2; 0 0];
phylogenyY = NaN(2,NsppY);
%phylogenyY(:,1) = [1,0];
phylogenyY(:,1:2) = [1 2; 0 0];

connectedness = zeros(1,T);
connectX = zeros(1,T);
connectY = zeros(1,T);

for n = 1:T % n is number of generations
    xall(n,:) = x(1,:);
    yall(n,:) = y(1,:);
    
    % we compute the average trait value of the interactors with any
    % particular species
    %avg value of the insects to which each fig is connected
    AvgXconnection = A* (y(1,iYextant))'./(sum(A'))';
    muI = AvgXconnection;
    %avg value of the figs to which each insect is connected
    AvgYconnection = A'*(x(1,iXextant))'./(sum(A))'; 
    muJ = AvgYconnection;
    Aold = A;
    % think of the X species as fig, basal spp
    
    
    wij = A;
    wji = A;
    w = A;
    wbari = ones(NextantsppX,1);
    wbarj = ones(NextantsppY,1);
    for i = 1:NextantsppX
        nI = sum(A(i,:));
        for j = 1:NextantsppY
            mJ = sum(A(:,j));
            wij(i,j) = ((theta+v)/(mJ*sqrt(theta^2+theta*v+v^2)))*exp(-1*(v*(x(1,iXextant(i))-muJ(j))*(x(1,iXextant(i))+muJ(j)-2*muI(i))+theta*(muI(i)-muJ(j))*(muI(i)+muJ(j)-2*x(1,iXextant(i))))/(theta^2+theta*v+v^2));
            wji(i,j) = ((theta+v)/(nI*sqrt(theta^2+theta*v+v^2)))*exp(-1*(v*(y(1,iYextant(j))-muI(i))*(y(1,iYextant(j))+muI(i)-2*muJ(j))+theta*(muJ(j)-muI(i))*(muI(i)+muJ(j)-2*y(1,iYextant(j))))/(theta^2+theta*v+v^2));
            w(i,j) = wji(i,j);
            %wibarj(i,j) = (1/mJ)*exp(-(x(1,iXextant(i))-muJ(j))*(x(1,iXextant(i))*v-2*muI(i)*(theta+v)+muJ(j)*(2*theta+v))/(2*(theta+v)^2));
            %wjbari(i,j) = (1/nI)*exp(-(y(1,iYextant(j))-muI(i))*(y(1,iYextant(j))*v-2*muJ(j)*(theta+v)+muI(i)*(2*theta+v))/(2*(theta+v)^2));
        end
    end
    
    %wi = sum(wij,2);%note wi(xi) = sum(wij(i,:))
    %wj = sum(wji,1);
    %wibar = sum(wibarj,2);
    %wjbar = sum(wjbari,1);
    %for i = 1:NextantsppX
    %    for j = 1:NextantsppY
    %       w(i,j) = wji(i,j)/(mJ*wjbar(j)); %relative fitness of insect i for a particular plant j  
   %     end
   % end
   
    m = 0;
    %Xspeciation = [];
    NoldsppX = NsppX;
    for i = iXextant %only look at the extant species
        % i is the index of each of the extant of individuals/species
        % recall we need to put NaN where spec have gone extinct, not
        % delete them completely, so this limits us to only calling extant
        % species to speciate and evolve
        m = m + 1; %counts which species we are on to use as an index for iXextant, A
                    %where i acts as an index for x, i.e iXextant(m) = i
      
        % variation from one generation to the next
        %
        nI = sum(Aold(m,:));
        thetaI = thetay*nI; %actually theta_I^2, 
        thetaI = thetay/nI;%CLT
        alpha = (thetaI^2+thetaI*v+v^2)/(thetaI^2+thetaI*v+v^2+thetax*v);
        gamma1 = thetaI*thetax/(thetaI^2+thetaI*v+v^2+thetax*v);
        x(2,i) = (alpha)*x(1,i)+(1-alpha)*AvgXconnection(m)+(AvgXconnection(m)-mean(AvgYconnection(find(Aold(m,:)))))*gamma1+randn*xix;
      
        
        % speciation event - which occurs with probability pX        
        if (rand < pX)
            NsppX = NsppX + 1; % keeps track of how many spp there are now
            %Xspeciation = [Xspeciation, i];  %stores indices of new species
                     
            if (NsppX > maxNX)
                disp('error, the number of species type X is greater than the number allocated')
                toobig = toobig+1;
                break
            end
            
            % add new spp trait
            x(2,NsppX) = x(1,i) + randn*xis;
            %A = [A; zeros(length(A(m,:)))]; does not inherit parental connections for
            % fig
            A = [A; A(m,:)]; %inherits parental connections
            w = [w; w(m,:)]; %and those connections inherit the same fitness value, at least initially
                        
            %now we make sure that we know that this new species is a
            %sister species to species i
            boo = find(phylogenyX(1,:)==i);
            for j = NsppX-1:-1:boo(1)+1
                phylogenyX(:,j+1) = phylogenyX(:,j);
            end
            phylogenyX(:,boo(1)+1) = [NsppX,n];            
        end
        
    end

    
    % think of the Y species as Fig Wasp
    m = 0;
    NoldsppY = NsppY;
    for i = iYextant %only look at the extant species
        % i is the index of each of the extant of individuals/species
        % recall we need to put NaN where spec have gone extinct, not
        % delete them completely, so this limits us to only calling extant
        % species to speciate and evolve
        m = m + 1; %counts which species we are on
      
        % then comes variation from one generation to the next
        mJ = sum(Aold(:,m));
        %thetaJ = mJ*thetax; %actually theta_J^2
        thetaJ = thetax/mJ; %CLT
        beta = (thetaJ^2+thetaJ*v+v^2)/(thetaJ^2+thetaJ*v+v^2+thetay*v);
        gamma2 = thetaJ*thetay/(thetaJ^2+thetaJ*v+v^2+thetay*v);
        y(2,i) = (beta)*y(1,i)+(1-beta)*AvgYconnection(m)+(AvgYconnection(m)-mean(AvgXconnection(find(Aold(:,m)))))*gamma2+randn*xix;
        
        
        % first we ask if there is a speciation event - which occurs with
        % probability p
        if (rand < pY)
            NsppY = NsppY + 1; % keeps track of how many spp there are now         
            if (NsppY > maxNY)
                disp('error, the number of species type Y is greater than the number allocated')
                toobig = toobig+1;
                break
            end
            
            % concatonate matrix w/the phylogeny of new spp
            y(2,NsppY) = y(1,i)+ randn*xis;
            A = horzcat(A, A(:,m));
            w = horzcat(w, w(:,m));
                    
            %now we make sure that we know that this new species is a
            %sister species to species i
            boo = find(phylogenyY(1,:)==i);
            for j = NsppY-1:-1:boo(1)
                phylogenyY(:,j+1) = phylogenyY(:,j);
            end
            phylogenyY(:,boo(1)+1) = [NsppY,n];
            end
        end       

    
    %connections!!
    %_________________________________________
    % interaction matrix - these are in the order in which spp appeared, 
    % but not including extinct species  
    count = 1; im = 0; jm = 0;
    iXextantnew = [iXextant (NoldsppX+1):NsppX]; %give the linages of X that are extant
    iYextantnew = [iYextant (NoldsppY+1):NsppY]; %give the linages of Y that are extant
    Anew = zeros(length(iXextantnew),length(iYextantnew));
    %if mod(n,10) == 0 
    %    Bnew = zeros(maxNY,maxNY);
    %end
    
    for i = iXextantnew
        im = im + 1;
        for j = iYextantnew
            jm = jm + 1;
            matching = c*A(im,jm)+ (1-c)*w(im,jm);
            if matching > rand
                Anew(im,jm)=1;
                %if mod(n,10) == 0
                %    Bnew(i,j)=1;
                %end
                count = count + 1;
            end
        end
        jm = 0;
    end
    
   % if NoldsppX~=NsppX
   %     A;
   %     Anew;
   % end
   %  if NoldsppY~=NsppY
   %     A
   %     Anew
   % end
    
    % check to see if this generation, there have been any changes to B
    % if not, move on, if there are, we need to compare the old interaction
    % matrix with the new interaction matrix and figure out what type of
    % change has been made
 
    %if (mod(n,10) == 0 && norm(Bnew - B)~=0)
    %if (norm(Bnew - B)~=0)
    %    k = find(Bnew-B);
    %    for i1 = 1:NsppX
    %        for i2 = i1:NsppX
    %            for j1 = 1:NsppY
    %                for j2 = j1:NsppY
    %                    Cnew = Bnew([i1 i2],[j1 j2]);
    %                    C = B([i1 i2],[j1 j2]);
    %                    if (Cnew-C)~=0
                            %now break into cases - duplication, sorting,
                            %and some cospeciation need speciation events to
                            %occur, so first we just look for host
                            %switching & co-speciation from duplicated lines, 
                            %then we only look for the others if a
                            %speciation event has occured
    %                        if C == [1 0; 0 0]
    %                            if Cnew == [1 1; 0 0]
    %                            end
    %                        end                      
    %                    end
    %                end
    %            end                
    %        end
    %    end
    %elseif (NsppX-NoldsppX>0 || NsppY-NoldsppY>0)
    %    sorting = sorting + NsppX-NoldsppX + NsppY-NoldsppY;
            %there was a speciation event, but no change in interaction
            %matrix, so the species didn't establish a connection and died
    %end

    
    %if count == 1
    %    recovery = recovery +1;
    %    if recovery==5
    %        disp('everything is extinct')
    %        break
    %    end
    %end
    
    %we now want to eliminate all species that have no connections to
    %anyone else
    %To look for rows of 0 do a sum(A') then a find
    XextantA = find(sum(Anew,2)); %these are all the x species that have 
    % at least one connection maintained to someone else
    iXextant = iXextantnew(XextantA); %everyone with non-zero connection has been eliminated
    NextantsppX = length(iXextant);
    
    %To look for columns of 0, do a sum(A) then a find
    YextantA = find(sum(Anew,1)); % these are all the Y species that survive 
    % that have connections to at least one other species in the other
    % class
    iYextant = iYextantnew(YextantA);%everyone with non-zero connection has been eliminated
    NextantsppY = length(iYextant);
       
    
    if (NsppX==0 && NsppY==0)
        disp('extinction of both species')
        extinctions = extinctions +1;
        break
    elseif NsppX==0
        disp('extinction of species X')
        extinctions = extinctions +1;
        break
    elseif NsppY==0
        disp('extinction of species Y')
        extinctions = extinctions +1;
        break
    end   
   
    %updating interaction matrix to exclude extinct species - they will
    %revert to the default NaN in the next generation
    A = Anew(XextantA, YextantA);
    %if mod(n,10) == 0
    %    B = Bnew;
    %    NoldsppXB = NsppX;
    %    NoldsppYB = NsppY;
    %end
    
    %[n, NextantsppX, NextantsppY]
    
    x(1,:) = x(2,:);
    x(2,:) = NaN(1, maxNX);
    y(1,:) = y(2,:);
    y(2,:) = NaN(1, maxNY);
    
    connectedness(1,n) = (sum((sum(A,1)),2))/(NextantsppX*NextantsppY);
    connectX(1,n)= (mean(sum(A,2)))/NextantsppY;
    connectY(1,n)= (mean(sum(A,1)))/NextantsppX;
    
    
    %if n == 100 || n==200 || n== 300
    %    disp(n)
    %    disp(phylogenyX)
    %    disp(phylogenyY)
    %    disp(iXextant)
    %    disp(iYextant)
    %    disp(A)
    %end
 
  end


%For node and edges plot
%iXextantpy = iXextant;
%for j = 1:NextantsppX
%iXextantpy (j) = find( phylogenyX(1,:)==iXextant(j));
%end
%iYextantpy = iYextant;
%for j = 1:NextantsppY
%iYextantpy (j) = find( phylogenyY(1,:)==iYextant(j));
%end

D=zeros(sum(sum(A)),2);
%Dx is 1st column, Dy is 2nd column
k=0;
for j=1:NextantsppY
    n=sum(A(:,j));
    for i=1:n
        k=k+1;
        D(k,2)=1/n;
    end
end

k=0;
for j=1:NextantsppX
    n=sum(A(j,:));
    for i=1:n
        k=k+1;
        D(j,1)=1/n;
    end
end

subplot(3,2,[2 4])
%plot(iXextantpy,ones(NextantsppX),'ro','MarkerSize',12,'MarkerFaceColor','r')
%hold on
%plot(iYextantpy,zeros(NextantsppY),'bo','MarkerSize',12,'MarkerFaceColor','b')
%ylim([-.5 1.5]);
%axis off
%maxextant=max([NsppY NsppX]);
%xlim([-.5 maxextant+.5]);
%hold on
    
%for i=1:NextantsppX
%    for j=1:NextantsppY
%        if A(i,j)==1
%            plot([iXextantpy(i),iYextantpy(j)],[1,0],'k')
%            hold on
%        end
%    end
%end
%hold off

iXextantpy = iXextant;
for j = 1:NextantsppX
iXextantpy (j) = find( phylogenyX(1,:)==iXextant(j));
end

iYextantpy = iYextant;
for j = 1:NextantsppY
iYextantpy (j) = find( phylogenyY(1,:)==iYextant(j));
end


Apy = zeros(size(A));
Axpy = (sortrows([1:NextantsppX; iXextantpy]',2))';
Aypy = (sortrows([1:NextantsppY; iYextantpy]',2))';
for i=1:NextantsppX
    for j=1:NextantsppY
        Apy(i,j)=1-A(Axpy(1,i),Aypy(1,j));
    end
end
imagesc(Apy)
%colormap(colorGray)
xlabel('Pollinator species')
ylabel('Plant species')

%plot each spp this generation
subplot(3,2,1)
plot([1:T],xall,'.', 'MarkerSize',5)
ylabel('Mean Trait x')
subplot(3,2,3)
plot([1:T],yall,'.', 'MarkerSize',5)
ylabel('Mean Trait y')
subplot(3,2,5)
%plot([1:T], connectY(1,:), [1:T], connectedness(1,:), [1:T], connectX(1,:))
plot([1:T], connectedness(1,:))
ylabel('Connectedness')
xlabel('Iterations')
subplot(3,2,6)
%hist(sum(A),0:max(sum(A)))
%Dependence
%hist(D)
%xlabel('Dependence')

%Connectivity distribution
maxC= max(sum(A,1));
hist(sum(A,1),1:maxC)
xlim([0,maxC])
xlabel('Number of connections per pollinator')

%visits = w;
%for i=1:size(w,1)
%    for j = 1:size(w,2)
%        visits(i,j) = round(sum(A(:,j),1)*w(i,j)+.40);
%    end
%end
%maxV=max(sum(visits,1));
%hist(sum(visits,1),0:maxV)
%xlabel('Number of plants visited per pollinator')

hold off
hold off
hold off


% Display counts for total number of speciation events
% later can input this into a file to be read later
disp('number of total speciation events:')
disp([NsppX-1 NsppY-1])
disp('number of extant species')
disp([NextantsppX NextantsppY])
disp('connectedness and average dependence')
disp([connectedness(1,T), mean(D(:,1)), mean(D(:,2))])

%disp('number of total sorting events:') 
%disp(sorting)

%print -djpeg90 -r0 CN_run16.jpg
%writematrix(Apy)
