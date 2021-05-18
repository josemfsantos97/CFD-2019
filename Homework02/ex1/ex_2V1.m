%% EX 2.1

close all
clear all
clc

%% Advection equation u_t + au_x = 0

c=0.1;                %Courant Number                
h=0.01;  %Grid Size
N=(2/h)+1;% Number of cells
x=-1:h:1;             %Centroid location
v_c=1;                %velocidade de conveção
dt=c*h/v_c;           %
j=0;         %critério de paragem

disp('Escolha o scheme:')
disp('1 - Explicit Second order Upwind (LUD)')
disp('2 - Implicit Second order Upwind (LUD)')
disp('3 - Explicit Minmod Flux Limiter')
disp('4 - Explicit QUICK')
disp('5 - Implicit QUICK')
disp('6 - Explicit First Order Upwind')
disp('7 - Implicit First Order Upwind')
disp('8 - Runge-Kutta 4th Order with upwind')
disp('9 - Runge-Kutta 4th Order with central')
metodo=input('');

%% Condição Inicial
u0=zeros(length(x),1);

a=0.5;
z=-0.7;
delta=0.005;
alpha=10;
beta=log(2)/(36*delta^2);

G1 = @(x) exp(-beta.*(x-z+delta).^2);
G2 = @(x) exp(-beta.*(x-z-delta).^2);
G3 = @(x) exp(-beta.*(x-z).^2);
H =  @(x) 1 - abs(10.*(x-0.1));
F1 = @(x) sqrt(1-(alpha.^2).*(x-a+delta).^2);
F2 = @(x) sqrt(1-(alpha.^2).*(x-a-delta).^2);
F3 = @(x) sqrt(1-(alpha.^2).*(x-a).^2);

%100.5 equivale a 0
u0((-0.8*(N/2)+(N/2)):(-0.6*(N/2)+(N/2)))= (1/6)*G1(-0.8:h:-0.6)+(1/6)*G2(-0.8:h:-0.6)+4/6*G3(-0.8:h:-0.6);
u0((-0.4*(N/2)+(N/2)):(-0.2*(N/2)+(N/2)))= 1;
u0((0*(N/2)+(N/2)):(0.2*(N/2)+(N/2)))= H(0:h:0.2);
u0((0.4*(N/2)+(N/2)):(0.6*(N/2)+(N/2)))=(1/6)*F1(0.4:h:0.6)+(1/6)*F2(0.4:h:0.6)+4/6*F3(0.4:h:0.6);

%% Explicit euler second order upwind
A = spdiags([(1-1.5*c)*ones(N,1),(2*c)*ones(N,1),(-0.5*c)*ones(N,1)],[0,-1,-2],N,N);
A(1:2,(N-1):N) = [-0.5*c ,(2*c);0 , -0.5*c]; % primeira e segunda linha dependem dos ultimos tmb

%% Implicit Euler second order Upwind
A1 = spdiags([(1+c+c/2)*ones(N,1),(-c-c)*ones(N,1),(c/2)*ones(N,1)],[0,-1,-2],N,N);
A1(1:2,N-1:N) = [c/2,(-c-c);0, c/2]; % depende do ultimo, periodic boundaries

%% Implicit Euler first order upwind
A2 = spdiags([(1+c)*ones(N,1),(-c)*ones(N,1)],[0,-1],N,N);
A2(1,N) = -c; % depende do ultimo, periodic boundaries

%% Explicit euler QUICK
A3 = spdiags([(-(3/8)*c)*ones(N,1),(1-(3/8)*c)*ones(N,1),((7/8)*c)*ones(N,1),(-(1/8)*c)*ones(N,1)],[1,0,-1,-2],N,N);
A3(1:2,(N-1):N) = [(-(1/8)*c) ,(7/8)*c ;0 , (-(1/8)*c)]; % primeira e segunda linha dependem dos ultimos tmb
A3(201,1)=(-(3/8)*c); %ultima linha depende do primeiro

%% Implicit Euler QUICK
A4 = spdiags([((3/8)*c)*ones(N,1),(1+(3/8)*c)*ones(N,1),(-(7/8)*c)*ones(N,1),((1/8)*c)*ones(N,1)],[1,0,-1,-2],N,N);
A4(1:2,(N-1):N) = [((1/8)*c) ,(-7/8)*c ;0 , ((1/8)*c)]; % primeira e segunda linha dependem dos ultimos tmb
A4(201,1)=((3/8)*c); %ultima linha depende do primeiro

u = real(u0);
x0=x;
phi = zeros(length(x),1);

switch metodo
    
    case 1
    %% Euler explicit second order     
while j<6
    phi=A*u;
    u=phi;
    x0=x0+dt;
    %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
    %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
    
end

    case 2
    %% Euler implicit second order
while j<6   
    
    phi=A1\u;
    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
    %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

    case 3
    %% Explicit Euler with MIN_MOD
while j<6
          
    for i=1:N
        %r(i)
        if i>1 && i<N
            ri=(u(i)-u(i-1))/(u(i+1)-u(i));
        elseif i==1
            ri=(u(1)-u(N))/(u(2)-u(1));
        elseif i==N 
            ri=(u(N)-u(N-1))/(u(1)-u(N));
        end
           
        %r(i-1)
        if i>2
            rii=(u(i-1)-u(i-2))/(u(i)-u(i-1));
        elseif i==2
            rii=(u(1)-u(N))/(u(2)-u(1));
        elseif i==1
            rii=(u(N)-u(N-1))/(u(1)-u(N));
        end
    
        psi_i=psi_minmod(ri);
        psi_ii=psi_minmod(rii);
        
        if i>1 && i<N
            phi(i) = u(i)*(1-c) + u(i-1)*(c) - FL(psi_i,psi_ii,u,i,N)*c;       
        elseif i==1
            phi(1) = u(1)*(1-c) + u(N)*(c) - FL(psi_i,psi_ii,u,i,N)*c;
        elseif i==N 
            phi(N) = u(N)*(1-c) + u(N-1)*c -FL(psi_i,psi_ii,u,i,N)*c;
        end  
        
    end

    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
    %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

    case 4
    %% Explicit Euler with QUICK
while j<6
   
    phi=A3*u;
    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
   %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

    case 5
    %% Implicit Euler with QUICK
while j<6   
    
    phi=A4\u;
    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
    %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

    case 6
    %% Explicit Euler with Upwind   
while j<6
    
    %para Euler explícito upwind => psi=0 e caso C=1 não há difusão, pelo
    %que não há erro de difusão
    
    for i=1:N
        
        if i>1 && i<=N
            phi(i) = u(i)*(1-c) + u(i-1)*(c);    
        elseif i==1
            phi(1) = u(1)*(1-c) + u(N)*(c);          
        end  
        
    end

    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
   %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

    case 7
    %% Implicit Euler with Upwind
while j<6   
    
    phi=A2\u;
    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
    %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end        

    case 8
%% Runge-Kutta 4th order with Upwind
while j<6
    u2=zeros(length(x),1);
    u3=zeros(length(x),1);
    u4=zeros(length(x),1);
    for i=1:N
        
        if i>1 && i<=N
            u2(i) = u(i) + u(i)*(-c*0.5) + u(i-1)*(c*0.5);    
        elseif i==1
            u2(1) = u(1) + u(1)*(-c*0.5) + u(N)*(c*0.5);          
        end  
    end
    
    for i=1:N
        if i>1 && i<=N
            u3(i) = u(i) + u2(i)*(-c*0.5) + u2(i-1)*(c*0.5);    
        elseif i==1
            u3(1) = u(1) + u2(1)*(-c*0.5) + u2(N)*(c*0.5);          
        end 
    end
    
    for i=1:N
        if i>1 && i<=N
            u4(i) = u(i) + u3(i)*(-c) + u3(i-1)*c;    
        elseif i==1
            u4(1) = u(i) + u3(1)*(-c) + u3(N)*c;          
        end
    end
    
    for i=1:N
        if i>1 && i<=N
            phi(i) = u(i) - (c/6)*( (u(i) - u(i-1)) + 2*( u2(i) - u2(i-1) ) + 2 * (u3(i) - u3(i-1)) + (u4(i) -u4(i-1)) );    
        elseif i==1
            phi(1) = u(1) - (c/6)*( u(1) - u(N) + 2*( u2(1) - u2(N) ) + 2 * (u3(1) - u3(N)) + u4(1) -u4(N) );         
        end 
    end

    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
   %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

    case 9
%% Runge-Kutta 4th order with Central
while j<6
    u2=zeros(length(x),1);
    u3=zeros(length(x),1);
    u4=zeros(length(x),1);
    
    
    
    for i=1:N
        
        if i>1 && i<N
            u2(i) = u(i) + (-c/2)*( ( u(i+1) - u(i-1) )/2 );    
        elseif i==1
            u2(1) = u(1) + (-c/2)*( ( u(2) - u(N) )/2 ); 
        elseif i==N
            u2(N) = u(N) + (-c/2)*( ( u(1) - u(N-1) )/2 ); 
        end  
    end    
    for i=1:N
        if i>1 && i<N
            u3(i) = u(i) + (-c/2)*( ( u2(i+1) - u2(i-1) )/2 );    
        elseif i==1
            u3(1) = u(1) + (-c/2)*( ( u2(2) - u2(N) )/2 );
        elseif i==N
            u3(N) = u(N) + (-c/2)*( ( u2(1) - u2(N-1) )/2 );
        end 
    end
    
    for i=1:N
        if i>1 && i<N
            u4(i) = u(i) + (-c)*( ( u3(i+1) - u3(i-1) )/2 );    
        elseif i==1
            u4(1) = u(1) + (-c)*( ( u3(2) - u3(N) )/2 );
        elseif i==N
            u4(N) = u(N) + (-c)*( ( u3(1) - u3(N-1) )/2 );
        end
    end
    
    for i=1:N   
        if i>1 && i<N
            phi(i) = u(i) - (c/6)*( ( ( u(i+1) - u(i-1) )/2 ) + 2*( ( u2(i+1) - u2(i-1) )/2 ) + 2 * ( ( u3(i+1) - u3(i-1) )/2 ) + ( ( u4(i+1) - u4(i-1) )/2 ) );    
        elseif i==1
            phi(1) = u(1) - (c/6)*( ( ( u(2) - u(N) )/2 ) + 2*( ( u2(2) - u2(N) )/2 ) + 2 * ( ( u3(2) - u3(N) )/2 ) + ( ( u4(2) - u4(N) )/2 ) );
        elseif i==N
            phi(N) = u(N) - (c/6)*( ( ( u(1) - u(N-1) )/2 ) + 2*( ( u2(1) - u2(N-1) )/2 ) + 2 * ( ( u3(1) - u3(N-1) )/2 ) + ( ( u4(1) - u4(N-1) )/2 ) );
        end 
    end

    u=phi;
    x0=x0+dt;
     %% Plot analytical solution
    if x0(1)>1.01
       x0=x;
       j=j+1;
    end
   %% Plot 1st revolution
    if j==1 && x0(1)==-1
        figure(1)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('1st Revolution')
        hold off
    end
    %% Plot 5th revolution
    if j==5 && x0(1)==-1
        figure(2)
        plot(x,real(u),'k')
        hold on
        plot(x0,real(u0),'r');
        title ('5th Revolution')
        hold off
    end
end

end

%%plot Error
figure()
u_exact=[real(u0(x0>1));real(u0(x0<=1))];
Erro=abs(u_exact-u);
plot(x,Erro)
title('Erro vs x')
xlabel('x')
ylabel('|Erro de difusão|')