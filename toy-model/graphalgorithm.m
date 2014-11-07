%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%                                                                                               %               
%                  Final Code - Algorithm Documentation                %
%                                                                                               %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   m-file name: graphalgorithm.m
%
%   code developper: Vassilis A. Stavrakas   
%   
%    e-mail: bill.scs7@gmail.com
%   
%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%                 Introduction              %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   This code implements a graph based approach to construct
%   a final signaling topology, combining experimental 
%   phosphoproteomic data with a PKN (Prior Knowledge 
%   Network).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Clearig the Workspace and the Command Window 
%   from previous work
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
%%
%%%%%%%%%%%%%%%%
%                               %
%   Data Reading       %
%                               %
%%%%%%%%%%%%%%%%
%
%   To read our network data we implement an algorithmic 
%   approach rather than ready MatLab functions.  We intent 
%   to make it easier for the user, because some reading functions 
%   as for example 'xlsread' do not always cooperate with all 
%   the Operational Systems (i.e Linux). This is a universal way 
%   to read every possible data format. 
%       
%%%%%%%%%%%%%%%
%   Network Data
%%%%%%%%%%%%%%%%%

netdata = importdata('Network_Generic.txt');

net = str2mat(netdata);
columns = 2;    % The possible columns number defined by the user

net2 = {length(net),columns};
row = 1;

for i = 1 : length(net)
    
    wordstart = 1;
    wordend = j;
    col = 1;
    
    net2{row, col} = net(i,1);
    
    for j = 1 : size(net,2)
                
        if ( strcmp(net(i,j),'	') == 0 )
            
           wordend  = j;
        
        else
            
            for  k = wordstart + 1 : wordend
              
                net2{row, col} = strcat(net2{row,col},net(i,k));
                
            end
            
            wordstart = j + 1;
            col = col + 1;
            net2{row, col} = net(i,j + 1);
            
        end
        
    end
    
    for k = wordstart + 1 : wordend
              
           net2{row, col} = strcat(net2{row,col},net(i,k));
                
    end
    
    row = row + 1;

end

network(:,1) = net2(:,1);
network(:,2) = net2(:,columns);

%%%%%%%%%%%%%%%%%%%%%%%%%
%   Stimuli - Signals Data Reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimuli = importdata('stimuli.txt');

signals = importdata('signal.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Reading the Experimental Matrix of 0 - 1
%   0: signifies unchanged state
%   1: signifies protein up regulation
% -1: signifies protein down regulation   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dependmtr = importdata('phosphos.txt');

% Making the dependencies matrix

dependencies={};
k=1;

for i=1:size(dependmtr,1)
 for j=1:size(dependmtr,2)
   
  if dependmtr(i,j)~=0
      
   dependencies{k,1}=stimuli{i}; 
   dependencies{k,2}=signals{j};
   k=k+1;
   
  end
  
 end
end

% Finding the unique sources and target nodes in the network

uniquesource=unique(dependencies(:,1));
uniquetarget=unique(dependencies(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clearing the unnecessary variables from the heap 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ('col','i','j','k','net','netdata', 'row','columns', 'wordend','wordstart');
clear ('textdata','net2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  %
%                       Counting the reactions                       % 
%                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   We count the frequence and we sort the reactors and we make the 
%   array of reactors incidence.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Reactors
%%%%%%%%%%%%%

reactors=network(:,1);  
reactors=sort(reactors);

match1=reactors{1};
react={ };
count1=1; 
k=1; 

for i=2:length(reactors)                      
  
  B=reactors{i};                                
  matching = strcmp(match1, B);                 
                                                 
  if matching==1                                 
     
      count1=count1+1;                            
  
  else
      
     react{k,1}=match1;
     react{k,2}=count1;                         
     match1=B;                                  
     k=k+1;                                     
     count1=1;                                  
  
  end                                           
                                                 
end    

%%%%%%%%%%%%%%%%%%%%%%%%                         
%          Last Array Element         %                         
%%%%%%%%%%%%%%%%%%%%%%%%

react{k,1}=match1;
react{k,2}=count1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Produced
%%%%%%%%%%%%%

produced=network(:,2);                          
produced=sort(produced);                         
                                                
match2=produced{1};                              
prod={ };                                       
count2=1;                                        
k=1;                                            
                                                 
for i=2:length(produced)      
    
  B=produced{i};                                
  matching = strcmp(match2, B);                 
                                                 
  if matching==1                                
      
      count2=count2+1;                           
  else
      
     prod{k,1}=match2;                          
     prod{k,2}=count2;                           
     match2=B;                                   
     k=k+1;                                     
     count2=1;                                  
     
  end                                          
                                                 
end  

%%%%%%%%%%%%%%%%%%%%%%%%                         
%          Last Array Element         %                         
%%%%%%%%%%%%%%%%%%%%%%%%

prod{k,1}=match2;                             
prod{k,2}=count2;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                               %     
%                  Finding Edge - Vertex number              %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Total Number of Nodes     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes=union(react(:,1),prod(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
%                         Graph Statistical Analysis                      %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                      

clc;

 fprintf(' ----------------------------------------\n');
fprintf('|                                                                         |\n');
fprintf('|                Graph Statistical Analysis                  |\n');
fprintf('|                                                                         |\n');
 fprintf(' ----------------------------------------\n');
 
fprintf('\n');
disp(['the total number of reactions is:' num2str(length(reactors))]);
fprintf('\n');
disp(['the total number of reactors is:' num2str(length(react))]);
fprintf('\n');
disp(['the total number of produced is:' num2str(length(prod))]);
fprintf('\n');
disp(['the total number of nodes in graph is:' num2str(length(nodes))]);
fprintf('\n');   
disp(['the total number of signals in graph is:' num2str(length(signals))]);
fprintf('\n');   
disp(['the total number of stimuli in graph is:' num2str(length(stimuli))]);
fprintf('\n'); 
disp(['the total number of source nodes in graph is:' num2str(length(uniquesource))]);
fprintf('\n'); 
disp(['the total number of targets in graph is:' num2str(length(uniquetarget))]);
fprintf('\n'); 

fprintf(' ----------------------------------------------------\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          
%                                       Adjacency   Matrix
%                                          Representation                                                      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   As it is known, in order to process graphs with a computer 
%   program, we first need to decide how to represent them within the 
%   computer. The first step in repersenting a graph is to map the 
%   vertex names to integers between 1 and  V. The main reason for 
%   this, is to make it possible to quickly access information 
%   corresponding to each vertex, using array indexing. Such a straight
%   forward representation for graphs is the so - called adjacency 
%   matrix representation. A V -by- V array of boolean values is 
%   maintained, a[i][j] set to 1, if there is an edge from vertex i to 
%   vertex j and 0 otherwise.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Making the Adjacency Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

adjmtrx=zeros(length(nodes));

 for i=1:length(network)

   r=network(i,1);
   p=network(i,2);
   
   for j=1:length(nodes)
       
    matching1=strcmp(nodes,r);
    matching2=strcmp(nodes,p);
    
   end
   
   pos1=find(matching1==1);
   pos2=find(matching2==1);
   
   adjmtrx(pos1,pos2)=1;
   
 end

% Adjmtrx verification 
%%%%%%%%%%%%%%%%%%%%%%%%%

% dokim=strcmp(nodes,'a20');
% ps1=find(dokim==1);
% dokim=strcmp(nodes,'traf6');
% ps2=find(dokim==1);
% adjmtrx(ps1,ps2)

%Delete the unecessary variables
clear('ans','ps1','ps2','dokim','matching1','matching2','pos1','pos2','r');
clear('k','match1','match2','count1','count2','B','p','i','j','ps','matching');
clear('textdata3');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                           Elementary Graph Algorithms                                 %         
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%           Strongly Connected Components       :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Finds the strongly connected components of the graph represented 
%   by matrix G using Tarjan's algorithm. 
%   A strongly connected component is a maximal group of nodes that 
%   are mutually reachable without violating the edge directions. 
%   Input G is an N-by-N sparse matrix that represents a graph. 
%   Nonzero entries in matrix G indicate the presence of an edge
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G=sparse(adjmtrx);

disp('the number of components found in the graph is:')
[S, C] = graphconncomp(G);  
disp(S)
disp('and to which component each node belongs')
disp(C);

%%%%%%%%%%%%%%%%%    
%       Spanning tree     :
%%%%%%%%%%%%%%%%%
%
%   A spanning tree must touch all the nodes and must be acyclic. 
%   G is an N-by-N sparse matrix whose lower triangle represents an 
%   undirected  graph. Nonzero entries in matrix G indicate the 
%   presence of an edge.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TF = graphisspantree(G);   %  returns logical 1 (true) if G is a spanning 
                                          %  tree, and logical 0 (false) 
    
if TF==1
    disp('the graph is spanning tree')
else
    disp('the graph is not a spanning tree')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum Spanning Tree :
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Finds an acyclic subset of edges that connects all the nodes in the 
%   undirected graph G and for which the total weight is minimized.
%   Weights of the edges are all n%the N-by-N sparse matrix G. 
%   Output Tree is a spanning tree represented by a sparse matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Tree, pred] = graphminspantree(G);

clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Transitive closure NOT Warshall's algorithm :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 

path = adjmtrx;
ind = length(adjmtrx);
counter=0;

 for i=1:ind
  for j=1:ind
   if path(i,j)==1
       for k=1:ind
         if path(j,k)==1
            path(i,k) = 1;
         end
       end
   end
  end
  
 end 
 
 % Nodes communication 

tran_clos={};
k=1;

for i=1:length(adjmtrx)
    
    tran_clos{k,1}=nodes{i};
    col=2;
    
    for j=1:length(adjmtrx)
        if (i~=j) && (path(i,j)==1)
            tran_clos{k,col}=nodes{j};
            col=col+1;
        end
    end
    
    k=k+1;
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%          Floyd - Warshall Algorithm   with Path Reconstruction            %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

react_new={};
prod_new={};
db_new=[];

for i=1:length(network)

  react_new{i}=network{i,1};
  prod_new{i}=network{i,2};
  db_new(i)= 1;

end

react_new=react_new';
prod_new=prod_new';
db_new=db_new';


    nodes_new=union(react_new,prod_new);

    dist =inf(length(nodes_new));
    next=zeros(length(nodes_new));
    B_new=nodes_new';
    
    for i=1:length(nodes_new)
        dist(i,i)=0;
    end
   
    for j=1:length(react_new)
        
      matching1=strcmp(B_new,react_new{j});  
      row_w=find(matching1==1);
      matching2=strcmp(B_new,prod_new{j}); 
      column_w=find(matching2==1);
      dist(row_w,column_w)=db_new(j);
      
    end

    for k=1:length(nodes_new)
      
        for i=1:length(nodes_new)
            for j=1:length(nodes_new)
                if dist(i,k)+dist(k,j)<dist(i,j)
                    dist(i,j)=dist(i,k)+dist(k,j);
                    next(i,j)=k;
                end
            end
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%
%       Path Reconstruction      %   
%%%%%%%%%%%%%%%%%%%%%%

   path_w=cell(length(nodes_new));

    for i=1:length(nodes_new)
        
      for j=1:length(nodes_new)
        path_w{i,j}=Path4(i,j,next);
      end
      
    end

 %Delete the unecessary variables
 clear ('col','counter','i','j','k','C','G','Tree','TF','S','ans','pred');
 clear('produced','reactors','source','target');
 
% Warshall's verification 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

source=dependencies(:,1);
target=dependencies(:,2);


for i=1:length(dependencies)
   
   A=source(i);
   B=target(i);
   dokim=strcmp(nodes_new,A);
   ps1=find(dokim==1);
   dokim=strcmp(nodes_new,B);
   ps2=find(dokim==1);
   d=dist(ps1,ps2);
   disp(d);   %  If dist between two nodes=1 --> there is no in-between  
                  %  path. There is only the reaction between them!!!
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          % 
%                           Checking the dependencies file                              %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Although the proteins connectivity in the final solution will be  a 
%  subset of the connectivity in the PKN,  we have to checkif the 
%  experimental data dictates connectivity that is not supported by the 
%  PKN.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(' ---------------------------------------------------\n');
fprintf('|                                                                                           |\n');
fprintf('|               Checking the dependencies file                            |\n');
fprintf('|                                                                                           |\n');
fprintf(' --------------------------------------------------\n');
 
 objective={};
 ct=0;
 
 for i=1:length(dependencies)
     
     A=source(i);
     B=target(i);
     
     matching=strcmp(A,nodes_new);
     pos1=find(matching==1);
     
     matching=strcmp(B,nodes_new);
     pos2=find(matching==1);
    
     if isempty(pos1)==1  
         
         fprintf('----------------------------------------------\n');
         fprintf('\n');
         disp(['The source node' cellstr(A)]);
         disp('does NOT exist in  the network data');
         fprintf('\n');
         fprintf('----------------------------------------------\n');
         ct=ct+1;
         objective{i}={};
         
     elseif isempty(pos2)==1  
         
         fprintf('----------------------------------------------\n');
         fprintf('\n');
         disp(['The target  node' cellstr(B)]);
         disp('does NOT exist in  the network data');
         fprintf('\n');
         fprintf('----------------------------------------------\n');
         ct=ct+1;
         objective{i}={};
         
     elseif isempty(pos1)==1 && isempty(pos2)==1
         
         fprintf('----------------------------------------------\n');
         fprintf('\n');
         disp([' The source node ' cellstr(A)]);
         disp([' and the target node ' cellstr(B)]);
         disp(' do NOT exist in the network data');
         fprintf('\n');
         fprintf('----------------------------------------------\n');
         ct=ct+1;
         objective{i}={};       
         
     else   
         
         objective{i}=path_w{pos1,pos2};
         
     end
    
 end
 
 
 objective=objective';
 
 if ct~=0
     
  disp([num2str(ct) ' dependencies do NOT exist in the network data ']);
 
 elseif ct==length(dependencies)
 
     disp( ' None dependency relation exists in the network data ');
 
 else
     
     disp( ' All the dependencies exist in the network data ');

end
%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                                                                                                         % 
 %                                       Direct Pathways                                          %
 %                                                                                                         %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % 
 %  We have already found the Direct Pathways.  As a  Direct Path 
 %  we call the path that results from Warshall's Algorithm. After  
 %  this, we check if this path is acceptable by seeing if each path-node 
 %  is a signal node. We find the non-empty paths positions in nepp 
 %  matrix and we check  if it is acceptable path.
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 fprintf(' --------------------------------------------------\n');
fprintf('|                                                                                          |\n');
fprintf('|                                    Direct Pathways                             |\n');
fprintf('|                                                                                          |\n');
 fprintf(' --------------------------------------------------\n');

 nepp=[];       % Array with non-empty-path-positions 
                     %  (direct path positions)
               
 rpp=[];         % Array with rest-path-positions
 k=1;
 m=1;
 
 for i=1:length(objective)
     
   if isempty(objective{i})==0
    
     nepp(k)=i;
     k=k+1;
    
   else
     
     rpp(m)=i;
     m=m+1;
     
   end
   
 end
   
 nepp=nepp';
 rpp=rpp';
 
                                       %%%%%%%%%%%%%%%% 
                                       %      Signal Nodes     %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  Here we find the index position of each signal node in the 
 %  nodes_new matrix             
 %
 
  signod=zeros(length(signals),1);
  
  for i=1:length(signals)
      
    A=signals(i);
    matching=strcmp(nodes_new,A);
    pos=find(matching==1);
    
    if isempty(pos)==1
        
    signod(i)=0;
    
    else
        
    signod(i)=pos;
    
    end
    
  end   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   %%%%%%%%%%%%%%%%%
                                   %     Stimuli Nodes      %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  Here we find the index position of each stimuli node in the 
 %  nodes_new matrix             
 %

  stimnod=zeros(length(stimuli),1);
  
  for i=1:length(stimuli)
      
    A=stimuli(i);
    matching=strcmp(nodes_new,A);
    pos=find(matching==1);
    
    if isempty(pos)==1
        
    stimnod(i)=0;
    
    else
        
    stimnod(i)=pos;
    
    end
    
  end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %           Presenting the Direct Paths               %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 % Here we present the Direct Warshall's Path  
 %
 
 finalgr={};
 fgct=1;
 ctr=0;
 
 for i=1:length(nepp)
   
   flag=0;  
   A=source(nepp(i));
   B=target(nepp(i));
   
   helpy=nodes_new(objective{nepp(i)});
   
   for j=1:length(helpy)
       
     compr2=strcmp(helpy{j,1},signals);
     pst2=find(compr2);
     compr1=strcmp(A,stimuli);
     pst1=find(compr1);
     
     if sum(compr2)==1 && dependmtr(pst1,pst2)==0
         flag=1;
     end
     
   end
   
   if flag==0
       
     disp([' From source ' cellstr(A)]);
     disp([' we re going to target ' cellstr(B) ]);
     disp('through the path: ');   
     fprintf('\n');
     disp (nodes_new(objective{nepp(i)}));
     fprintf('------------------------------------------------\n');
     fprintf('\n');
     
     helpa={};
     bl=1;
     helpa(bl)=A;
     
     for k=1:length(helpy)
      
         bl=bl+1;     
         helpa{bl}=helpy{k};
     
     end
     
     helpa(bl+1)=B;
     helpa=helpa';
     
     
     for k=1:length(helpa)-1
        
        finalgr{fgct,1}=helpa{k};
        finalgr{fgct,2}=helpa{k+1};
        finalgr{fgct,3}='objective';
        fgct=fgct+1;
        
     end
     
   else
       
     ctr=ctr+1;       
     rpp(length(rpp)+1)=nepp(i);
   
   end
   
 end
 
 %Delete the unecessary variables
 clear('A','B','matching','pos1','pos2','i','ps1','ps2','dokim','d','match');
 clear('data','db_new','nodes','path','prod','react','textdata','tran_clos');
 clear('ans','k','j');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                                      No Possible Pathways                                  %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%  Here we find the NO possible pathways for the dependencies we 
%  have set. As no possible pathway we call  the path through 
%  which we will never reach the target node because it aint exist. And 
%  if there is NO Warshall's path, we seek other possible paths through 
%  the source nodes. If there is NO other possible path, then there is 
%  NO path at all.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(' --------------------------------------------------\n');
fprintf('|                                                                                          |\n');
fprintf('|                               No possible Pathways                         |\n');
fprintf('|                                                                                          |\n');
 fprintf(' -------------------------------------------------\n');


rlbmtr=zeros(length(rpp),1);  % Reliability Matrix shows if there is a 
                                              %  possible path --> 1  or  NOT --> 0 !!!

for i=1:length(rpp);
 
  s=source(rpp(i));
  t=target(rpp(i));
 
  mtcs=strcmp(s,nodes_new);
  nds=find(mtcs==1);
 
  mtct=strcmp(t,nodes_new);
  ndt=find(mtct==1);
 
  d=dist(nds,:);
  
  if length(unique(d))>2
      
     rlbmtr(i,1)=1;
  
  end
 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Presenting the NO possible paths               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Here we present the NO paths because there is no Warshall's path 
%  and not any other possible path through other nodes
%

 for i=1:length(rlbmtr)
     
   if rlbmtr(i)==0
    
    disp([' From source ' cellstr(source(rpp(i)))]);    
    disp('there is NO posible path: ');  
    fprintf('\n');
    disp(['to target ' cellstr(target(rpp(i))) ]);    
    disp(' in this network !!! ');
    fprintf('\n');
    fprintf('-------------------------------------------------\n');
    fprintf('\n'); 

   end
   
 end
%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %                                                                                                         %
 %                                           Alternative Paths                                    %
 %                                                                                                         %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %  So far, we have already found the Direct Warshall Paths without 
 %  an intermediate signal node, and we have also found the NO 
 %  possible paths from each source node to the respective target  
 %  (based on the dependencies). So, here we examine the rest 
 %  dependencies:
 %
 %    a) If there is a direct path (NOT the shortest one) through some 
 %        other nodes then we have a desirable path. If this path goes 
 %        through a signal node then is not permitable and we face what 
 %        we call "CONFLICT", because we can ONLY go to the target node 
 %        through another signal node.
 %
 %    b) If there is no other direct pathm then we face the similar 
 %        case of NO possible pathways.
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 fprintf(' --------------------------------------------------\n');
fprintf('|                                                                                          |\n');
fprintf('|                                  Alternative Paths                              |\n');
fprintf('|                                                                                          |\n');
 fprintf(' --------------------------------------------------\n');

%   Finding the array with the possible alternative pathways positions 
%   for the rest source nodes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=1; 
pos=[];

for i=1:length(rpp);
  
  if rlbmtr(i)==1
      
    s=source(rpp(i));
    t=target(rpp(i));
 
    mtcs=strcmp(s,nodes_new);
    nds=find(mtcs==1);
 
    mtct=strcmp(t,nodes_new);
    ndt=find(mtct==1);
 
    d=dist(nds,:);
    uni=unique(d);
    
    l=1;
    for j=1:length(uni)
     
     if uni(j)~=0 && uni(j)~=inf   
        
       for k=1:length(d)
           
         if dist(nds,k)==uni(j)
             
           pos(m,l)=k;    %   Rows are the rest nodes - columns are the 
                                 %    alternative dist positions
           
           l=l+1;            % Next column
           
         end
         
       end
       
     end
     
    end

  
  m=m+1;                % Next row
  end 
 
end

%
%    Finding the  indexes of the dependencies nodes     
%                      in the nodes_new matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  sourcenod=zeros(length(source),1);
  targetnod=zeros(length(target),1);
  
  for i=1:length(dependencies)
       
      A=source(i);
      match1=strcmp(A,nodes_new);
      
      if sum(match1)~=0
          
          sourcenod(i)=find(match1==1);
      
      end
      
      B=target(i);
      match2=strcmp(B,nodes_new);
      
      if sum(match2)~=0
      
          targetnod(i)=find(match2==1);
      
      end
      
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the rest array with the source nodes with existing alter path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 rest=[];
 j=1;
 for i=1:length(rlbmtr)
   
   if rlbmtr(i)==1
       
     rest(j,1)=rpp(i);
     j=j+1;
     
   end
   
 end
 
% Presenting the alternative paths for the rest source nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 rest1=rest;
 for i=1:size(pos,1)
      
  for j=1:nnz(pos(i,:))   
   
   if targetnod(rest(i))~=0  
    pth=path_w{pos(i,j),targetnod(rest(i))};
    
    ct=0;
    if isempty(pth)==0
 
      for k=1:size(pth,2)

         if sum(pth(1,k)==signod)==0 
             
           ct=ct+1;
           
         else
             
           posit2=find(pth(1,k)==signod);
           a=(sourcenod(rest(i))==stimnod);
           posit1=find(a);
           
           if dependmtr(posit1,posit2) ~= 0 
               
             ct=ct+1;
                
           end
           
         end
         
      end
    end

   end
   
  if sourcenod(rest(i))~=0
      
   pth1=path_w{sourcenod(rest(i)),pos(i,j)};
    
    ct1=0;
    if isempty(pth1)==0
        
      for k=1:size(pth1,2)

         if sum(pth1(1,k)==signod)==0
             
           ct1=ct1+1;
           
         else
             
           posit2=find(pth1(1,k)==signod);
           a=(sourcenod(rest(i))==stimnod);
           posit1=find(a);
           
           if dependmtr(posit1,posit2) ~=0
               
             ct1=ct1+1;
                
           end
        
         end
         
      end
      
    end
  end
   
   if isempty(pth1)==0 &&  isempty(pth)==0
       rest1(i)=inf;    
    if ct==size(pth,2) && ct1==size(pth1,2)
     
       fprintf('\n');
       disp([' From source: ' cellstr(source(rest(i)))]);
       disp([' we go to target: ' cellstr(target(rest(i)))]);
       disp(' through the path:' );
       disp(nodes_new(pth1));
       disp(nodes_new(pos(i,j)));
       disp(nodes_new(pth));
       fprintf('\n');
       fprintf('-----------------------------------------------\n');
       rest1(i)=0;
       
     helpa={};
     bl=1;
     helpa(bl)=cellstr(source(rest(i)));
     
     for k=1:length(nodes_new(pth1))
      bl=bl+1;
      helpa(bl)=nodes_new(pth1(k));
     end
     
     bl=bl+1;
     helpa(bl)=nodes_new(pos(i,j));
     
     for k=1:length(nodes_new(pth))
      bl=bl+1;
      helpa(bl)=nodes_new(pth(k));
     end
     
     helpa(bl+1)=cellstr(target(rest(i)));
     helpa=helpa';
     
     for k=1:length(helpa)-1
        
        finalgr{fgct,1}=helpa{k};
        finalgr{fgct,2}=helpa{k+1};
        finalgr{fgct,3}='alternative';
        fgct=fgct+1;
        
     end

    end
    
   end
      
  end
   
 end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                                           Finding Conflicts                                     %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Here we find the conflicts. If there is a path with a signal node
%  through it, then we have a conflict. If there is NO possible path 
%  neither an alternative path, then we present the NO existing paths. 
%  We check if the Path_w of rest1 is empty. If it is true, then there is 
%  NO possible pathway.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
fprintf(' --------------------------------------------------\n');
fprintf('|                                                                                          |\n');
fprintf('|                                    Finding Conflicts                            |\n');
fprintf('|                                                                                          |\n');
fprintf(' --------------------------------------------------\n');
 
 rest2=rest1;
 for i=1:length(rest1)
     
   ctcf=0;  
   ctcf1=0;
   ctcf2=0;
   
  if rest1(i)~=0
      
   if sourcenod(rest(i))~=0 && targetnod(rest(i))~=0
        
    helpyp=path_w{sourcenod(rest(i)),targetnod(rest(i))};
    helpy=nodes_new(helpyp);
    
    if isempty(helpyp)==1 
        
     disp([' From source ' cellstr(source(rpp(i)))]);    
     disp('there is NO posible path: ');  
     fprintf('\n');
     disp(['to target ' cellstr(target(rpp(i))) ]);    
     disp(' in this network !!! ');
     fprintf('\n');
     fprintf('------------------------------------------------\n');
     fprintf('\n'); 
     rest2(i)=0;
     
    else 
     
     for k=1:size(helpyp,2)
      
      posit2=find(helpyp(k)==signod);
      a=(sourcenod(rest(i))==stimnod);
      posit1=find(a);   
           
      if sum(helpyp(k)==signod)==1 && dependmtr(posit1,posit2)==0 
         
       disp([' From source: ' cellstr(source(rest(i)))]);
       disp([' we go to target: ' cellstr(target(rest(i)))]);
       disp(' through the path:' );
       disp(nodes_new(helpyp));
       fprintf('\n'); 
       disp('but the signal conflict is '); 
       disp(nodes_new(helpyp(k)))   
       fprintf('\n');
       fprintf('-----------------------------------------------\n');
       rest2(i)=1;
       
       if (sum(sum(strcmp(nodes_new(helpyp(k)),finalgr)))==0)
            ctcf=ctcf+1;
       end
       
      end
      
     end
     
    end
    
   end 
     
  elseif rest1(i)==inf
   
   if sourcenod(rest(i))~=0
       
    pth1=path_w{sourcenod(rest(i)),pos(i,j)}; 
    
    posit2=find(pth1(k)==signod);
    a=(sourcenod(rest(i))==stimmnod);
    posit1=find(a);   
    
    for k=1:size(pth1,2)

      if sum(pth1(k)==signod)==1 && dependmtr(posit1,posit2)==0   
         
       disp([' From source: ' cellstr(source(rest(i)))]);
       disp([' we go to intermediate node: ' cellstr(pos(i,j))]);
       disp(' through the path:' );
       disp(nodes_new(pth1));
       disp(nodes_new(pth));
       fprintf('\n'); 
       disp('but the signal conflict is '); 
       disp(nodes_new(pth1(k)))   
       fprintf('\n');
       fprintf('-----------------------------------------------\n');
       rest2(i)=1;
       
       if (sum(sum(strcmp(nodes_new(pth1(k)),finalgr)))==0)
            ctcf1=ctcf1+1;
       end
       
      end
      
    end
    
   end
   
   if targetnod(rest(i))~=0
       
    pth=path_w{pos(i,j),targetnod(rest(i))}; 
    
    
    for k=1:size(pth,2)
        
    posit2=find(pth(k)==signod);
    a=(sourcenod(rest(i))==stimmnod);
    posit1=find(a);  

      if sum(pth(k)==signod)==1 && dependmtr(posit1,posit2)==0   
        
       disp([' From intermediate node: ' cellstr(pos(i,j))]);
       disp([' we go to target: ' cellstr(target(rest(i)))]);
       disp(' through the path:' );
       disp(nodes_new(pth1));
       disp(nodes_new(pth));
       fprintf('\n'); 
       disp('but the signal conflict is '); 
       disp(nodes_new(pth(k)));   
       fprintf('\n'); 
       fprintf('-----------------------------------------------\n');
       rest2(i)=1;

       if (sum(sum(strcmp(nodes_new(pth(k)),finalgr)))==0)
            ctcf2=ctcf2+1;
       end
       
      end
      
    end
    
   end
  
  end
  
  if ctcf >= 0
       
     helpa={};
     bl=1;
     helpa(bl)=source(rest(i));
     
     for k=1:length(helpyp)
      bl=bl+1;
      helpa{bl}=helpy{k};
     end
     
     helpa(bl+1)=target(rest(i));
     helpa=helpa';
     
     
     for k=1:length(helpa)-1
        
        finalgr{fgct,1}=helpa{k};
        finalgr{fgct,2}=helpa{k+1};
        finalgr{fgct,3}='conflict';
        fgct=fgct+1;
        
     end
     
  end
  
    if ctcf1+ctcf2 >= 0
      
     helpa={};
     bl=1;
     helpa(bl)=cellstr(source(rest(i)));
     
     for k=1:length(nodes_new(pth1))
      bl=bl+1;
      helpa(bl)=nodes_new(pth1(k));
     end
     
     for k=1:length(nodes_new(pth))
      bl=bl+1;
      helpa(bl)=nodes_new(pth(k));
     end
     
     helpa(bl+1)=cellstr(target(rest(i)));
     helpa=helpa';
     
     
     for k=1:length(helpa)-1
        
        finalgr{fgct,1}=helpa{k};
        finalgr{fgct,2}=helpa{k+1};
        fgct=fgct+1;
        
     end
     
    end
    
 end

%Delete the unecessary variables
 clear('A','B','B_new','ans','column_w','ct','ct1','ctr','d','flag','helpy');
 clear('helpyp','i','ind','j','k','l','m','match1','match2','row_w','rest');
 clear('s','t','uni','pth','pth1','pth2','matching1','matching2','mtcs');
 clear('mtct','nds','ndt')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                                    Presenting the Conflicts                                 %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Here, we present finally the conflict nodes in the pathways.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ct=0;
for i=1:length(rest2)
    
 if rest2(i)==1
     
   ct=ct+1;
   
 end
 
end

fprintf('\n');
disp([' " We have ' num2str(ct) ' Conflicts in the network " ']);
fprintf('\n');

%Delete the unecessary variables
clear('ct','i');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                                        Exporting the Reaults                                %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(finalgr)
  
 if isempty(finalgr{i,1})==0 && isempty(finalgr{i,2})==0
     
  A=finalgr{i,1};
  B=finalgr{i,2};
  
 for j=1:length(finalgr)
   if j~=i && isempty(finalgr{j,1})==0 && isempty(finalgr{j,2})==0
       
       if (strcmp(A,finalgr{j,1})==1 && strcmp(B,finalgr{j,2})==1)
  
         finalgr{j,1}={};
         finalgr{j,2}={};
         
       end
       
   end
 end
 
 end

end

fstimnodes=zeros(length(stimuli),1);

for i=1:length(stimuli)
  for j=1:length(finalgr)
      
    A=stimuli(i,1);
    B=finalgr(j,1);
    
    if strcmp(A,B)==1
        
        fstimnodes(i,1)=1;
        
    end
    
  end
end

% 
% Merging the finalgr array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnet2={};

k=1;

for i = 1:length(finalgr)
    
    if isempty(finalgr{i,1}) == 0 && isempty(finalgr{i,2}) == 0
    
        A = finalgr{i,1};
        B = finalgr{i,2};
        
    end
    
    for j = 1: length(network)
        
        C = network{j,1};
        D = network{j,2};
        
        if strcmp(A,C) == 1 &&  strcmp(B,D) == 1

            fnet2{k,1} = A;  
            fnet2{k,2} = B;  
            fnet2{k,3} = j;
            k = k+1;
            
        end
        
    end
    
end
    
        
% Reactions Clustering
 %%%%%%%%%%%%%%%%%%%%%%%
 
for i = 1 :length(fnet2)
   for j = i +1 : length(fnet2)
       
       if fnet2{i,3} == fnet2{j,3}
       
                fnet2{j,3} = [];
           
       end
       
   end
end
 
fnet = {};
k = 1;

for i = 1:length(fnet2)
    
    if isempty(fnet2{i,3}) == 0
        
            fnet{k,1} = fnet2{i,1};
            fnet{k,2} = fnet2{i,2};
            fnet{k,3} = fnet2{i,3};
            k = k +1;
            
    end
    
end
          
% File Opening
%%%%%%%%%%%%%%%%

p=fopen('fnet.dot','w');

if p== 0
  fprintf('Error in opening a file\n');
else
fprintf('Succesfully opened the file\n');    
end

fprintf(p,'digraph G{');
fprintf(p,'\n');

% Reactions Writing
%%%%%%%%%%%%%%%%%%%%

fprintf(p,'// Reactions (Edges)');
fprintf(p,'\n');  


 for j=1:length(network)
  
  flag=0;
 
  for i=1:length(fnet)
   
   if isempty(fnet{i,3})==0
       
    if fnet{i,3}==j    
   
     flag=1;
     pos=i;
     break
    
    end
    
   end
   
  end
 
  if flag==1
      
      fprintf(p,' "%s" -> "%s" [penwidth=2.00,color="black"]; \n',fnet{pos,1},fnet{pos,2}); 
  
 else    
     fprintf(p,' "%s" -> "%s" [style="dashed"]; \n',network{j,1},network{j,2});
     
 end
  
 end 
    
% Nodes Writing
%%%%%%%%%%%%%%%%%
 
fprintf(p,'// Compounds (Nodes)');
fprintf(p,'/n');  


nodes=nodes_new;

for i=1:length(nodes)
 
  A=nodes(i);  
  flag=0;
  
  for j=1:length(fnet) 
  
   C=fnet{j,1};   
   D=fnet{j,2};   
  
   if (strcmp(A,C)==1|| strcmp(A,D)==1)   
    
       flag=1;
       pos=i;
       break
        
   end
  
  end
     
  if flag==1
      
       fprintf(p,'"%s" [shape=oval,style="filled,rounded",color=black,fillcolor=cyan] \n', nodes{i,1});
              
  else
       
       fprintf(p,'"%s" [shape=oval,style="dashed"] \n', nodes{i,1});
  
   end
   
end

% Stimuli Writing
%%%%%%%%%%%%%%%%%%

stimuli = importdata('stimuli.txt');

for i=1:length(stimuli)
   
  flag=0;  
  A=stimuli(i,1);
  for j=1:length(fnet)
     
    if isempty(fnet{j,3})==0  
      
      if (strcmp(A,fnet{j,1})==1)
        fprintf(p,'"%s" [shape=box3d,color=dodgerblue4, style="filled,rounded", fillcolor=cornflowerblue] \n', stimuli{i,1});
        flag=1;
      end
      
    end
    
  end
  
  if flag==0
      
    fprintf(p,'"%s" [shape=box3d,style="filled,rounded",fillcolor=gray] \n', stimuli{i,1});  
      
  end
  
end
   
% Signals Writing
%%%%%%%%%%%%%%%%%%

signals = importdata('signal.txt');

for i=1:length(signals)
  
  flag=0;
  B=signals(i,1);
  for j=1:length(fnet)
    
    if isempty(fnet{j,2})==0 && isempty(fnet{j,3})==0
    
      if (strcmp(B,fnet{j,2})==1)
         fprintf(p,'"%s" [shape=oval,color=black, style="filled,rounded", fillcolor=chartreuse] \n', signals{i,1});
         flag=1;
      end
    
    end
    
  end
   
  if flag==0
      
    fprintf(p,'"%s" [shape=oval,style="filled,rounded",fillcolor=gray] \n', signals{i,1});  
   
  end
  
end

fprintf(p,'// Logical ANDs (Dots) \n');
fprintf(p,'}');
fclose all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                               Final Network Reactions Counter                        %                
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('\n');
fprintf(' --------------------------------------------------\n');
fprintf('\n');
fprintf('\n');
disp([' Reactions in the final network : ', num2str(length(fnet)) ]);
fprintf('\n');
fprintf('\n');
fprintf(' --------------------------------------------------\n');
fprintf('\n');
fprintf('\n');

%Delete the unecessary variables
 clear('A','B','C','ans','D','ct','a','ctr','bl','compr1','compr2');
 clear('p','posit1','posit2','pst1','pst2','flag','i','j','k','pos','rest1','rest2');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                      %
%                              Comparison  to  ILP method                            %                
%                                                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Reading ILP  solution network
%%%%%%%%%%%%%%%%%%%%%%%%%

netdata = importdata('ilp.txt');

net = str2mat(netdata);
columns = 2;    % The possible columns number defined by the user

net2 = {length(net),columns};
row = 1;

for i = 1 : length(net)
    
    wordstart = 1;
    wordend = j;
    col = 1;
    
    net2{row, col} = net(i,1);
    
    for j = 1 : size(net,2)
                
        if ( strcmp(net(i,j),'	') == 0 )
            
           wordend  = j;
        
        else
            
            for  k = wordstart + 1 : wordend
              
                net2{row, col} = strcat(net2{row,col},net(i,k));
                
            end
            
            wordstart = j + 1;
            col = col + 1;
            net2{row, col} = net(i,j + 1);
            
        end
        
    end
    
    for k = wordstart + 1 : wordend
              
           net2{row, col} = strcat(net2{row,col},net(i,k));
                
    end
    
    row = row + 1;

end

teonet(:,1) = net2(:,1);
teonet(:,2) = net2(:,2);
  
%  First Estimation
%%%%%%%%%%%%%%%%%%%
%
%  As a first comparison - estimation method of both graphs we can  
%  use reactions number's division in each final graph  
%
 if length(fnet) > length(teonet)
     
   estim = (length(teonet)/length(fnet))*100 ;
   
 else
     
   estim = (length(fnet)/length(teonet))*100 ;

 end

fprintf('\n');
fprintf('\n');

 disp([' A first similarity estimation between the two graphs is:  ' num2str(estim) '%']);

fprintf('\n');
fprintf('\n');


% Jaccard Index
%%%%%%%%%%%%%%%%%%%
%
% The Jaccard index, also known as the Jaccard similarity coefficient  
% is a statistic used for comparing the similarity and diversity of 
% sample sets. The Jaccard coefficient measures similarity between  
% sample sets, and is defined as the size of the intersection divided  
% by the size of the union of the sample sets.
%
%%%%%%%%%%%%%%%%%%%%
%                              δυο ψεματα
%    Graph Solution Matrix     
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
bill=[ ] ;
k=1;

for i=1:length(fnet)
    
    A=fnet{i,1};
    B=fnet{i,2};
    
    for j=1:length(network)
        
      if strcmp(A,network{j,1})==1 && strcmp(B,network{j,2})==1
          
          bill(k)=j ;
          k=k+1 ;
          
      end
      
    end
    
end
 
%%%%%%%%%%%%%%%%%%%%
%                               
%    ILP  Solution Matrix      
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

teodor=[ ] ;
k=1;

for i=1:length(teonet)
    
    A=teonet{i,1};
    B=teonet{i,2};
    
    for j=1:length(network)
        
      if strcmp(A,network{j,1})==1 && strcmp(B,network{j,2})==1
          
          teodor(k)=j ;
          k=k+1 ;
          
      end
      
    end
    
end

Jaccard = [length(intersect(bill,teodor)) / length(union(bill,teodor))]*100;


 fprintf(' ------------------------------------------------\n');
 fprintf('|                                                                                      |\n');
 fprintf('|                                     Jaccard Index                            |\n');
 fprintf('|                                                                                      |\n');     
 fprintf(' ------------------------------------------------\n');
 
 fprintf('\n');
fprintf('\n');
 disp([' The Jaccard similarity index between the two graphs is:  ' num2str(Jaccard) '%']);
fprintf('\n');
fprintf('\n');
%%

 fprintf(' --------------------------------------------------\n');
 fprintf('|                                                                                         |\n');
 fprintf('|                                         THE  END                                 |\n');
 fprintf('|                                                                                         |\n');     
 fprintf(' --------------------------------------------------\n');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                          %
%                                                   END                                                %
%                                                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
