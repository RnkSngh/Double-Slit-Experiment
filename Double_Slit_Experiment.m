%creates a video of a 2D simulation of the double slit experiment
%Numerically solves the u_dot_dot- c^2*grad(div(u))=f equation for u, 
%with given neumann boundary conditions of c^2*grad(u)*n=h on the neuman nodes,
%and dirichlet conditions of u=2*A*pi*omega*cos(2*pi*omega*t), and
%initial conditions of u, udot= 0 for t=0
%Uses triangular elements, and linear basis functions

clear

%USER DEFINED VARIABLES
timestep= .0001; %increasing timesteps cause longer runtimes
theta=0.5;
mass_lump= false; %change to true if you want to have mass lumping; only for theta=0
tf=3; %the total integration time
number_of_frames = 400; %this is to see how many frames we can use to input in the video ********** Requires the final time/timestep to be greater than number of trames
output_name= 'Theta=0.5'; %the output file that this will be named to

%GIVEN VARIABLES
A=.05; %amplitude of the dirichlet condition
omega=3; 
c=1;
h=0;
g=@(t) A*sin(2*pi*omega*t);

%first solve for the g deravitives, then input into function handles
g_dot= @(t) A*2*pi*omega*cos(2*pi*omega*t);
g_dot_dot= @(t) -4*A*(pi^2)* (omega^2)* sin(2*pi*omega*t);

frame_reset= ceil(tf/(timestep*number_of_frames)); %frame reset is the number of timesteps between recording each frame

%graph times are when a timestep should be graphed
times= 0:timestep:tf;

%Br is the gradients of basis functions in local coordinates
Br = [-1 1 0; -1 0 1];

%load given coordinates, connectivity, and boundary data
connectivity= load("dsg-connectivity.dat");
coordinates=  load("dsg-coordinates.dat");
dirichlet= load("dsg-dirichlet.dat");
neumann= load("dsg-neumann.dat");

%calculate number of nodes and number of elements
nn= length(coordinates(:, 1));
ne= length(connectivity(:, 1));

%output will keep track of what will be recorded
output= zeros(nn, number_of_frames);

%load video writer
v=VideoWriter(output_name , 'MPEG-4');
vFramerate= 15;
open(v);


%for solving at the non_dirichlet nodes, we keep track of non_dirichlet
%index
non_dirichlet= setdiff(1:nn, dirichlet);

%now we get an array of non_dirichlet nodes for both u and v
non_dirichlet_nodes=[non_dirichlet, non_dirichlet+nn];



%initialize K, M, F
K= sparse(nn, nn);
M= sparse(nn, nn);
F_last= zeros(nn, 1);
F_current= zeros(nn, 1);

%iterate through, adding each Ke, Me, and initial Fe
for element_index= 1:ne
    %get the coordinates for all nodes in the element
    nodeIDs=connectivity(element_index ,:);
    global_coordinates= coordinates(nodeIDs, :);
    
    %split into x and y coordinates
    x= global_coordinates(:,1);
    y= global_coordinates(:,2);
    
    %compute the jacobian, its inverse, and the determinite
    %of basis functions deravitvies wrt global coordinates
    J=[x(2)-x(1), x(3)-x(1); y(2)-y(1), y(3)-y(1)];
    invJ= inv(J);
    det_J= det(J);
    
    %add the Ke to each element; this is a constant so we take this out of
    %the integral and multiply by the area of a triangle (.5)
    K(nodeIDs, nodeIDs)= K(nodeIDs,nodeIDs)+.5*Br'*invJ*invJ'*Br*det_J;
    M(nodeIDs, nodeIDs)= M(nodeIDs, nodeIDs) + det_J*[2 1 1; 1 2 1; 1 1 2]/24;
    
end

%we store matricies called Ksum and Msum so that we can use in the
%definition of F when we are iterating. This is so that we can decrease the
%number of calculations which need to be done 
K_sum= zeros(nn, 1);
M_sum= zeros(nn, 1);

for i=nn
    K_sum(i)= sum(K(i, :));
    M_sum(i)= sum(M(i, :));
end


%if we are mass lumping, we sum all entries in m
if mass_lump==true
    for i=1:nn
        row_sum= sum(M(i, :));
        M(i, :)=0;
        M(i, i)= row_sum;
    end
end


%the Left hand side does not depend on time so we can iterate through
LHS= [M, -theta*timestep*M; theta*timestep*K,  M];

%calculate the part of right hand side which is constant wrt time
constant_RHS= [M, (1-theta)*timestep*M; (theta-1)*timestep*K, M];

%apply dirichlet condition and given initial conditions for u and v
u=zeros(2*nn, 1);
u(dirichlet)= g(0); 
u(dirichlet+nn)= g_dot(0); %u consists of [u;v], and its length is 2*nn
u(non_dirichlet_nodes)= 0;
d=u;

%add to output to graph
output(:, 1)= d(1:nn);

Nt=ceil(tf/timestep); % the number of timesteps that we will iterate throguh 


%now iterate through time steps to build ds
for time_index= 1:Nt %time index keeps track of which time we have
    time= times(time_index);
    

    d_last=d; %store d of the last value as d_last
    F_last=F_current; %store the last F value into a vector
    d=zeros(2*nn, 1); %reset the current d so that we can calculate it later
    

    g_value=g(time)*theta + (1-theta)*g(time-timestep); %uses the theta method to approximate between time steps; this is used for dirichlet conditions
    g_dot_value= g_dot(time)*theta + (1-theta)*g_dot(time-timestep); %this value will be used in the dirichlet conditions
    g_dot_dot_value= theta*g_dot_dot(time) + (1-theta)*g_dot_dot(time-timestep); %this value will be used for calculating F
    
    if tf==6 & time> 3 %checks for last part 
        g_value=0;
        g_dot_value=0;
        g_dot_dot_value=0;
    end
        
    %build F
    for row= 1:nn
        F_current(row)= -(K_sum(row)*g_value) -(M_sum(row)*g_dot_dot_value);
    end
    
    %build the RHS
    RHS= constant_RHS*d_last+ [zeros(nn, 1); (theta*timestep*F_current +((1-theta)*timestep*F_last))];
    
    d(dirichlet)= g_value; 
    d(dirichlet+nn)= g_dot_value; %d consists of [u;v], and its length is 2*nn
    
    d(non_dirichlet_nodes)= LHS(non_dirichlet_nodes, non_dirichlet_nodes)\RHS(non_dirichlet_nodes);
    
    %check if we need to graph, then add to output
    if mod(time_index, frame_reset)==0
        print_string = sprintf ('%d percent complete ' , round(time*100/tf)) ;
        output(:, time_index/frame_reset)= d(1:nn);
        disp ( print_string ) ;
    end
        
end

   
%% DISPLAY RESULTS
trisurf(connectivity,coordinates(:,1),coordinates(:,2),output(:,end))
l = [];
Frame(tf/timestep) = struct('cdata',[],'colormap',[]);
fig=figure;
for ti = 1:number_of_frames
    trisurf(connectivity,coordinates(:,1),coordinates(:,2),output(:,ti))
    view(20,40);
    zoom out
    axis([-.5, 1, 0, 1, -0.2, 0.2]);
    if ti*frame_reset>length(times)
        close(v);
        return
    end
    
    t=times(ti*(frame_reset));
    title("u(x,t) for theta=0.5, timestep=.001, t=" + string(times(ti*(frame_reset))))
    xlabel("x1");
    ylabel("x2");
    drawnow;
    Frame(ti) = getframe(fig);
    writeVideo(v, Frame(ti)); 
end

close(v);

