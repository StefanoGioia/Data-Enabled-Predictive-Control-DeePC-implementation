%% DeePC 2x2, 2 inps 2 outs %%
%Stefano Gioia, 2020 

%Problem: 2D speed tracking with front & side accelerations as inputs

clc
clear all
close all

global x 

%Initial state (front and side speed)
% 
x = 5;    
x(2,1) = 3.5; 

%% time delta
global dt
dt= 0.01; 

%scaling!

%% Reference trajectory mission

%Traj. generation
% [r1;r2] = [ sinusoidal speed on front, constant on side]

%Data to insert
steps= 66;
spd= 10; %to have an order of magnitude of speeds
nofwaves= 3;
beginWave = 0.3; % spd*sin(beginWave) = initial side speed
%

r= spd*ones(1,steps);
r(2,:) = spd*sin(beginWave: (nofwaves*2*pi)/(steps-1) :nofwaves*2*pi+beginWave);
%plot(r(2,:))

%% Set up the problem 

usize= 2;
ysize= 2;
% Tini >= l(lag)>= n(=#states) = here= 2 ,  Tfut= Tfuture= predicted steps

% min columns of Hankel matrix UorY : T-(Tini+Tfut)+1 >= #ctrls*(Tini+Tfut) 
% -> T>= (#ctrl+1)*(Tini+Tfut)-1 >= 14 (with Tini=2)
% I use 19 to be safe 

T = 19; %Higher can improve results
Tini = 3; 

% Tfuture = predicted steps
predicted_steps = 3 ;

applied_steps = 2;


hor_sideHankel=  T - (Tini+predicted_steps)+1 ;

ver_sideTotHankel = (usize+ysize)*(Tini+predicted_steps);   

%Persistently exciting condition applied on U

Hu = zeros(usize*(Tini+predicted_steps),hor_sideHankel);
 while rank(Hu) ~= usize*(Tini+predicted_steps)
     
     disp('Computing Hankel...')

% % Generate (initial) history - past data
%Attempting to make a persistently exciting input: row rank!
% u= 0.005*(rand(usize,T)-0.5).*( (10*ones(usize,T)) .^ (10*rand(usize,T)-5) ); 
    
    u =  900*(rand(usize,T)-0.5); %0.005*(rand(usize,T)-0.5).*( (10*ones(usize,T)) .^ (10*rand(usize,T)-5) ); %
    %900 multiplier found based on results: works better than with 1 
    
    %build output history
    for i = 1:T
        
        y(:,i) = SpeedModel(u(:,i)) ;
        
    end
    
    
    %build Hankels
    
    Hp = zeros( (usize+ysize)*Tini , hor_sideHankel);
    
    for j= 1:Tini %+1: first counts
        Hp(usize*(j-1)+1 :usize*j,:) = u(:, j : j+  hor_sideHankel-1) ;  % U % [U;Y] % (i-1): 1 timestep shift
          Up(usize*(j-1)+1 :usize*j,:) = u(:, j : j+  hor_sideHankel-1) ;
        Hp(Tini*usize+ ysize*(j-1)+1 :Tini*usize+ysize*j,:) = y(:, j: j+hor_sideHankel-1) ;  % Y % [U;Y]
        
    end
    
    Hf = zeros( (usize+ysize)*predicted_steps,hor_sideHankel);
    for j= 1: predicted_steps
        
        Hf(usize*(j-1)+1 :usize*j,:) = u(:,Tini+ j: Tini+ j+hor_sideHankel-1) ;  % U % [U;Y] % (i-1): 1 timestep shift
        Uf(usize*(j-1)+1 :usize*j,:) = u(:, Tini+ j: Tini+ j+hor_sideHankel-1) ;  % U % [U;Y] % (i-1): 1 timestep shift
        
        Hf(predicted_steps*usize+ ysize*(j-1)+1 :predicted_steps*usize+ ysize*j,:) = y(:, Tini+ j: Tini+ j+hor_sideHankel-1) ;  % Y % [U;Y]
        Yf(ysize*(j-1)+1 :ysize*j,:) = y(:, Tini+ j: Tini+ j+hor_sideHankel-1) ;  % Y % [U;Y]
    end
    
    Htot = [Hp; Hf];
    Hu = [Up;Uf];        
 end
 
 disp('Hankel rank ok')

%% Optimization set up

%Equality constr.
A = Hp ;

%Cost
Qcost= zeros(hor_sideHankel) ;
% if R=1; Q=1;R=2*R; Q=2Q; then
%u
Q = 0.003;
R = 1e-9*Q;
sqUf = sum(Uf.^2);
sqYf = sum(Yf.^2);
for k= 1: hor_sideHankel
    
    Qcost(k,k) = R*sqUf(k) + Q*sqYf(k);
    
end

%separate components
for j= 1: predicted_steps
    Ufcomp1(j,:) = Uf((j-1)*usize+1,:);
    Ufcomp2(j,:) = Uf((j-1)*usize+2,:);
end
for j= 1: predicted_steps
    Yfcomp1(j,:) = Yf((j-1)*ysize+1,:);
    Yfcomp2(j,:) = Yf((j-1)*ysize+2,:);
end

%Mixed products due to elevation
triang = zeros(hor_sideHankel) ;
for j= 1: (hor_sideHankel-1)
    
    for k = (j+1): hor_sideHankel
        
        for kk= 1:predicted_steps
            
            triang(j,k) = triang(j,k)+ R*Ufcomp1(kk,j)*Ufcomp1(kk,k) + R*Ufcomp2(kk,j)*Ufcomp2(kk,k) ;
            
            triang(j,k) = triang(j,k)+ Q*Yfcomp1(kk,j)*Yfcomp1(kk,k) + Q*Yfcomp2(kk,j)*Yfcomp2(kk,k) ;
            
            
        end
        triang(k,j) = triang(j,k);
    end
    
end

Qcost = 2*(Qcost+ triang); % 2*: quadprog divides by two


%% Predictive Control

gstory= zeros(hor_sideHankel,1);
fvalstory=nan;

cycles = floor((steps-(predicted_steps-applied_steps))/applied_steps);

for i = 1: cycles
    
      
    %solve for g
    
    %Build Hankel matrices from past data and last behavior prediction
   
        % p : past, lastStepPast = T-predicted_steps;
                
% % No need to update Hankel
        
        % Optimization set up
                
        %Equality constr.
        %b to update? Yes! A piece of trajectory that should have continuity
        
        b=zeros((Tini)*(usize+ysize),1);
        b(1:Tini*usize,1) = reshape(u(:,end-Tini+1:end),[Tini*usize,1]);
        b(Tini*usize+ 1:Tini*usize+ Tini*ysize,1) = reshape(y(:,end-Tini+1:end),[Tini*ysize,1]);
         %Use Tini as first part of behavior
%        b(1:Tini*usize,1) = reshape(u(:,1:Tini),[Tini*usize,1]);
%       b(Tini*usize+ 1:Tini*usize+ Tini*ysize,1) = reshape(y(:,1:Tini),[Tini*ysize,1]);
        
        %Linear cost        
        
        Lcost= zeros(hor_sideHankel,1);
        
%         sumyf1 =sum(Yfcomp1);
%         sumyf2 =sum(Yfcomp2);
        
        rYf1 = r(1,(i-1)*applied_steps+1:(i-1)*applied_steps+predicted_steps) * Yfcomp1; 
        rYf2 = r(2,(i-1)*applied_steps+1:(i-1)*applied_steps+predicted_steps) * Yfcomp2; 
        
        
        for j= 1:hor_sideHankel
            
            Lcost(j) = -2*Q*(rYf1(j) + rYf2(j));
                      
        end
                

%          qpoptions.ConstraintTolerance =1e-3;
%         qpoptions.MaxIterations = 600;
%         qpoptions.OptimalityTolerance =1e-5;
%          qpoptions.StepTolerance=1e-9;
        
   %Optimize
%     [g,fval] = quadprog(Qcost,Lcost,[],[],A,b,[],[],[],qpoptions); %solver stopped prematurely: cannot meet suddenly reference
       [g,fval] = quadprog(Qcost,Lcost,[],[],A,b);
       %store g also !
       
       gstory(:,i)  = g;
       fvalstory(i) = fval;
       
       newus= reshape(Uf*g, [usize, predicted_steps]) ;
       
       newys= reshape(Yf*g, [ysize, predicted_steps]) ;
       
       %Apply steps
       
       for ii = 1:applied_steps
           
           y(:,end+1) = SpeedModel(newus(:,ii)) ;
           
       end
       u(:,end+1:end+applied_steps) = newus(:,1:applied_steps) ; 
       
end

computed_steps = applied_steps*cycles;

%see results
figure
subplot(3,1,1)
plot(y(2,T+1:T+computed_steps),'r')
hold on
subplot(3,1,2)
plot(r(2,1:computed_steps),'k')
hold on
subplot(3,1,3)
plot(y(2,T+1:T+computed_steps)-r(2,1:computed_steps) )
%plot((y(2,T+1:T+computed_steps)-r(2,1:computed_steps))./r(2,1:computed_steps) )
hold off
%comp1
figure
subplot(3,1,1)
plot(y(1,T+1:T+computed_steps),'r')
hold on
subplot(3,1,2)
plot(r(1,1:computed_steps),'k')
hold on
subplot(3,1,3)
plot(y(1,T+1:T+computed_steps)-r(1,1:computed_steps) )




%% Physics model

%Inputs:front, side acceleration [u1;u2]--Outputs (as states): front, side speed [y1,y2]

function out = SpeedModel(U)
    
    global x dt
    out = x(:,end) + U*dt ;    
    x(:,end+1) = out ;      % x tracks out history
    
end

