%% DeePC 2x2, 2 inps 2 outs%%
%Stefano Gioia, 2020 %

%Problem: 2D speed tracking with front & side accelerations as inputs

clc
clear all
close all
% States
global x 

%Initial state (front and side speed)
x = 5;    
x(2,1) = 3.5;


%% time delta
global dt
dt= 0.001; 

%% Reference trajectory mission

%Traj. generation
% [r1;r2] = [ sinusoidal speed on front, constant on side]

%Data to insert
steps= 66;
spd= 10; %to have an order of magnitude of speeds
nofwaves= 3;
beginWave = 0.3; % spd*sin(beginWave) = initial side speed
%

global r

r= spd*ones(1,steps);
r(2,:) = spd*sin(beginWave: (nofwaves*2*pi)/(steps-1) :nofwaves*2*pi+beginWave);
%plot(r(2,:))

%% Set up the problem 

global usize ysize T predicted_steps hor_sideHankel

usize= 2;
ysize= 2;
% Tini >= l(lag)>= n(=#states) = here= 2 ,  Tfut= Tfuture= predicted steps

% Tfuture = predicted steps
predicted_steps = 1 ;

applied_steps = 1;

% min columns of Hankel matrix UorY : T-(Tini+Tfut)+1 >= #ctrls*(Tini+Tfut) 
% -> T>= (#ctrl+1)*(Tini+Tfut)-1 >= 14 (with Tini=2)
% I use 50 to be safe 
T = 50; 

Tini = 4 ; %here states=2, so put 4 to be safe

% sideHankel= (T+1)/2; 
hor_sideHankel= T - (Tini+predicted_steps)+1 ;

% % Generate (initial) history - past data
%Attempting to make something really persistently exciting, but fit to the
%magnitude of the reference trajectory : it's important!
%Rank-> check row rank


u= 0.005*(rand(usize,T)-0.5).*( (10*ones(usize,T)) .^ (10*rand(usize,T)-5) );

%build output history
for i = 1:T
    
    y(:,i) = SpeedModel(u(:,i)) ;
    
end

%% Predictive Control

global Uf Yf 
gstory= zeros(hor_sideHankel,1);
fvalstory=nan;

cycles = floor((steps-(predicted_steps-applied_steps))/applied_steps);

global i

for i = 1:cycles
    
    %I apply s=1 of the predicted steps
    
    %Reminder: in each iteration, there could be iterations to converge to solution
    % %solve for g
    
    
    if i==1  %Hankel(model) not to be changed
        
        %Build Hankel matrices from past data and last behavior prediction
        
        % p : past
        %lastStepPast = T-predicted_steps;
        
        Hp = zeros( (usize+ysize)*Tini , hor_sideHankel);
        
        for j= 1:Tini %+1: first counts
            Hp(usize*(j-1)+1 :usize*j,:) = u(:,applied_steps*(i-1)+ j  :applied_steps*(i-1)+ j+  hor_sideHankel-1) ;  % U % [U;Y] % (i-1): 1 timestep shift
            
            Hp(Tini*usize+ ysize*(j-1)+1 :Tini*usize+ysize*j,:) = y(:,applied_steps*(i-1)+ j:applied_steps*(i-1)+ j+hor_sideHankel-1) ;  % Y % [U;Y]
            
        end
        
        Hf = zeros( (usize+ysize)*predicted_steps,hor_sideHankel);
        for j= 1: predicted_steps
            
            Hf(usize*(j-1)+1 :usize*j,:) = u(:,applied_steps*(i-1)+Tini+ j:applied_steps*(i-1)+Tini+ j+hor_sideHankel-1) ;  % U % [U;Y] % (i-1): 1 timestep shift
            Uf(usize*(j-1)+1 :usize*j,:) = u(:,applied_steps*(i-1)+Tini+ j:applied_steps*(i-1)+Tini+ j+hor_sideHankel-1) ;  % U % [U;Y] % (i-1): 1 timestep shift
            
            Hf(predicted_steps*usize+ ysize*(j-1)+1 :predicted_steps*usize+ ysize*j,:) = y(:, applied_steps*(i-1)+Tini+ j: applied_steps*(i-1)+Tini+ j+hor_sideHankel-1) ;  % Y % [U;Y]
            Yf(ysize*(j-1)+1 :ysize*j,:)  = y(:, applied_steps*(i-1)+Tini+ j: applied_steps*(i-1)+Tini+ j+hor_sideHankel-1) ;  % Y % [U;Y]
            
            %Rank-> check row rank (assumed ok here)
            
        end
        
        Htot = [Hp; Hf];
        
        % Optimization set up
        
        %Equality constr.
        A = Hp;
    end
    %update for initial condition
        b=zeros((Tini)*(usize+ysize),1);
        b(1:Tini*usize,1) = reshape(u(:,end-Tini+1:end),[Tini*usize,1]);
        b(Tini*usize+ 1:Tini*usize+ Tini*ysize,1) = reshape(y(:,end-Tini+1:end),[Tini*ysize,1]);

        % to give dimensions to fmincon 
        g0 = rand(hor_sideHankel,1)-0.5*ones(hor_sideHankel,1);

        [g,fval] = fmincon(@DeePCcost,g0,[],[],A,b);

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


function cost = DeePCcost(xx)
      
    global r predicted_steps ysize Yf Uf i
    R=0.00000001; Q=1e11*R;
        
    %   Q and R here just multiplicator 1
    cost= Q*sum((Yf*xx - reshape(r(:,i:i+predicted_steps-1),[ysize*predicted_steps,1])).^2)+ R*sum((Uf*xx).^2);
     
    
end
