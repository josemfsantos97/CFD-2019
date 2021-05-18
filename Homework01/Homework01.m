clc
clear all

%% Problem Parameters
M=0.9;
airfoil_start=2;
c=1;
slope=0.06;

%% Number of cells
cellx=150;
celly=150;

%% Grid size
xdomain = 5;
ydomain = 5;

%% Problem Constants
Msub=0.8;
p_inf=10^5;
T_inf=288;
R=287;
rho_inf=p_inf/(R*T_inf);
gamma=1.4;
T=288;
R=287;
C=(1-M.^2);
Uinf=M*sqrt(gamma*T*R);

%% Cell size
dx=xdomain/cellx;
dy=ydomain/celly;

dxdy=dx/dy;
dydx=dy/dx;
%% Variables
x=cellx;
y=celly;
total=x*y;

if M<=0.8
    [sol,u,v,value]=subsonic_field(M,xdomain,ydomain,cellx,celly,airfoil_start,c);
elseif M>=1.2
    supersonic_field(M,xdomain,ydomain,cellx,celly,airfoil_start,c,slope);
elseif M>0.8 && M<0.95
    %% Get inital solution for phi and dphi/dx
    [sol,u,v,value,coef]=subsonic_field(Msub,xdomain,ydomain,cellx,celly,airfoil_start,c);

    %% Matrix A
    k=1:1:y-2;
    iter=0;
    error=1;
    tic
    while iter<150
        coef=sparse(total,total);
        iter=iter+1;
        display(num2str(iter))
        for i=1:1:total
            flag=0;
            K(i,1)=(1-M^2-(1+gamma)*(M^2/Uinf)*u(i));
            if i==1                 %%Bottom-Left
                if K(i)>0
    %                 display('bottom-left positivo',num2str(i))
                    coef(i,i)=-(3*K(i)*dydx+1*dxdy);
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i+x)=1*dxdy;
                    flag=1;
                else
    %                 display('bottom-left negativo',num2str(i))
                    coef(i,i)=(2*K(i)*dydx-dxdy);
                    coef(i,i+x)=dxdy;
                    flag=1;
                end
            end
            if i==2                                         %%Bottom-Left 2nd cell
                if K(i)>0
    %                 display('bottom-left 2 positivo',num2str(i))
                    coef(i,i)=-(2*K(i)*dydx+1*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i+x)=1*dxdy;
                    flag=1;
                else
    %                 display('bottom-left 2 negativo',num2str(i))
                    coef(i,i)=K(i)*dydx-dxdy;
                    coef(i,i-1)=-3*K(i)*dydx;
                    coef(i,i+x)=dxdy;
                    flag=1;
                end
            end

            if i>2 && i<x                                   %%Bottom cells
                if K(i)>0               
    %                 display('bottom positivo',num2str(i))
                    coef(i,i)=-(2*K(i)*dydx+1*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i+x)=1*dxdy;   
                    flag=1;
                else 
    %                 display('bottom negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-dxdy);
                    coef(i,i-1)=-2*K(i)*dydx;
                    coef(i,i-2)=K(i)*dydx;
                    coef(i,i+x)=dxdy;
                    flag=1;
                end 
            end
            if i==x                                         %Bottom-Right
                if K(i)>0
    %                 display('bottom-right positivo',num2str(i))
                    coef(i,i)=-(3*K(i)*dydx+1*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+x)=1*dxdy;
                    flag=1;
                else
    %                                 display('bottom right negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-dxdy);
                    coef(i,i-1)=-2*K(i)*dydx;
                    coef(i,i-2)=K(i)*dydx;
                    coef(i,i+x)=dxdy;
                    flag=1;
                end
            end
            if ismember(i,1+k*x)== 1                %Middle-Left Cells
                if K(i)>0
    %                 display('middle-left positivo',num2str(i))
                    coef(i,i)=-(3*K(i)*dydx+2*dxdy);
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i+x)=1*dxdy;
                    coef(i,i-x)=1*dxdy;
                    flag=1;
                else
    %                 display('middle-left negativo',num2str(i))
                    coef(i,i)=(2*K(i)*dydx-2*dxdy);
                    coef(i,i+x)=dxdy;
                    coef(i,i-x)=dxdy;
                    flag=1;
                end
            end
            if ismember(i,2+k*x)==1                    %Middle-Left 2nd Cells
                if K(i)>0 
    %                 display('middle-left 2 positivo',num2str(i))
            coef(i,i)=-(2*K(i)*dydx+2*dxdy);
            coef(i,i-1)=K(i)*dydx;
            coef(i,i+1)=K(i)*dydx;
            coef(i,i-x)=1*dxdy;
            coef(i,i+x)=1*dxdy;
                    flag=1;
                else
    %                 display('middle-left 2 negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-2*dxdy);
                    coef(i,i-1)=-3*K(i)*dydx;
                    coef(i,i+x)=dxdy;
                    coef(i,i-x)=dxdy;
                    flag=1;
                end
            end
            if ismember(i,x+k*x)==1                     %Middle-Right Cells
                if K(i)>0
    %                 display('middle-right positivo',num2str(i))
                    coef(i,i)=-(3*K(i)*dydx+2*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+x)=1*dxdy;
                    coef(i,i-x)=1*dxdy;
                    flag=1;
                else
    %                 display('middle-rght negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-2*dxdy);
                    coef(i,i-1)=-2*K(i)*dydx;
                    coef(i,i-2)=K(i)*dydx;
                    coef(i,i-x)=dxdy;
                    coef(i,i+x)=dxdy;
                    flag=1;
                end
            end
            if i==1+x*(y-1)                     %%Top-Left Cell
                if K(i)>0
    %                 display('top-left positivo',num2str(i))
                    coef(i,i)=-(3*K(i)*dydx+1*dxdy);
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i-x)=1*dxdy;
                    flag=1;
                else
    %                 display('top-left nagativo',num2str(i))
                    coef(i,i)=(2*K(i)*dydx-dxdy);
                    coef(i,i-x)=dxdy;
                    flag=1;
                end
            end
            if i==2+x*(y-1)                     %%Top-Left 2nd Cell
                if K(i)>0
    %                 display('top-left 2 positivo',num2str(i))
                    coef(i,i)=-(2*K(i)*dydx+1*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i-x)=1*dxdy;
                    flag=1;
                else
    %                 display('top-left 2 negativo',num2str(i))
                    coef(i,i)=K(i)*dydx-dxdy;
                    coef(i,i-1)=-3*K(i)*dydx;
                    coef(i,i-x)=dxdy;
                    flag=1;
                end
            end

            if i>2+x*(y-1) && i<total                 %%Top Cells
                if K(i)>0
    %                 display('top positivo',num2str(i))
                    coef(i,i)=-(2*K(i)*dydx+1*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i-x)=1*dxdy;
                    flag=1;
                else
    %                 display('top negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-dxdy);
                    coef(i,i-1)=-2*K(i)*dydx;
                    coef(i,i-2)=K(i)*dydx;
                    coef(i,i-x)=dxdy;
                    flag=1;
                end
            end
            if i==total                         %%Top-Right
                if K(i)>0
    %                 display('top right positivo',num2str(i))
                    coef(i,i)=-(3*K(i)*dydx+1*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i-x)=1*dxdy;
                    flag=1;
                else
    %                 display('top right negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-dxdy);
                    coef(i,i-1)=-2*K(i)*dydx;
                    coef(i,i-2)=K(i)*dydx;
                    coef(i,i-x)=dxdy;
                    flag=1;
                end
            end
            if flag==0                      %%Middle Cells
                if K(i)>0
    %                 display('middle positivo',num2str(i))
                    coef(i,i)=-(2*K(i)*dydx+2*dxdy);
                    coef(i,i-1)=K(i)*dydx;
                    coef(i,i+1)=K(i)*dydx;
                    coef(i,i-x)=1*dxdy;
                    coef(i,i+x)=1*dxdy;
                else
    %                 display('middle negativo',num2str(i))
                    coef(i,i)=(K(i)*dydx-2*dxdy);
                    coef(i,i-1)=-2*K(i)*dydx;
                    coef(i,i-2)=K(i)*dydx;
                    coef(i,i-x)=dxdy;
                    coef(i,i+x)=dxdy;
                end
            end

        end
        if iter==1
                sol_old=sol;
        end

%             sol_new = qmr(coef,value,1e-8,2000); 
            sol_new = coef\value;
            solution=0.1*sol_new+0.9*sol_old;
%             solution=sol_new;
            error=coef*solution-value;
            error_max(iter)=max(error);
            u=u_field(solution,x,y,dx,K);
            display(num2str(error_max));
            sol_old=solution;
    end
    toc
    solution=vec2mat(solution,x);
    iter_vector=[1:1:iter];
    %% Plot of Potential Field
    figure (4)
    X=linspace(dx,xdomain-dx,cellx);
    Y=linspace(dy,ydomain-dy,celly);
    surf(X,Y,solution,'EdgeColor','none')
    view(0,90)
    colorbar
    
    %% Calculations
    Utot =u+Uinf;                       %Total Velocity X
    Velocity=sqrt(Utot.^2+v.^2);        %Velocity Magnitude
    T0=T_inf*(1+(gamma-1)/2*M^2);       %Stagnation Temperature

    %Temperature Field
    T_field=(T0*R*gamma/(gamma-1)-0.5*Velocity.^2)*(gamma-1)/(gamma*R);

    %Mach Number Field
    Mach_local=Velocity./((gamma*R*T_field).^0.5);

    %Pressure Field
    P_field=p_inf*(T_field/T_inf).^(gamma/(gamma-1));

    %Stagnation Pressure Field
    P0_field=P_field.*(1+((gamma-1)/2)*Mach_local.^2).^(gamma/(gamma-1));

    %Density Field
    rho_field=rho_inf*(P_field/p_inf).^(1/gamma);
    
    
    %% Plot of Mach Flow Field
    figure(5)
    surf(X,Y,vec2mat(Mach_local,x),'EdgeColor','none');
    view(0,90);
    colorbar
    
    figure(6)
    plot(iter_vector,error_max)
    
else
    display('Impossible to solve for Mach between 0.96 and 1.2')
end
function [v] = v_field(sol,x,y,dy,dx,Uinf,airfoil_start,c)
    total=x*y;
    v=zeros(total,1); %central differences
    for i=1:1:total
        if i<=x
            if dx/2+(i-1)*dx>=airfoil_start && dx/2+(i-1)*dx<airfoil_start+c/2 %Airfol ascent
                v(i)=0.06*Uinf;
            elseif dx/2+(i-1)*dx>airfoil_start+c/2 && dx/2+(i-1)*dx<=airfoil_start+c %Airfoil descent
                v(i)=-0.06*Uinf;
            else
                v(i)=(sol(i+x)-sol(i))/dy;
            end
        elseif i>=1+x*(y-1)
            v(i)=(sol(i)-sol(i-x))/dy;
        
        elseif i>x && i<1+x*(y-1)
            v(i)=(sol(i+x)-sol(i-x))/(2*dy);
        end
    end
end
function [u] = u_field(sol,x,y,dx,K)
    %% Compute U
    total=x*y;
    u=zeros(total,1); %central differences
    k=0:1:y-1;
    for i=1:1:total
        if ismember(i,1+k*x)==1
            u(i)= (sol(i))/dx;
        elseif ismember(i,x+k*x)==1
            u(i)= (sol(i)-sol(i-1))/dx;
        else       %Subsonico
            u(i)= (sol(i+1)-sol(i-1))/(2*dx);
%         else                %SuperSonico
%             u(i)= (sol(i)-sol(i-1))/dx;
        end
    end
end

function supersonic_field(M,xdomain,ydomain,cellx,celly,airfoil_start,c,slope)
    %% Problem Constants
    p_inf=10^5;
    T_inf=288;
    R=287;
    rho_inf=p_inf/(R*T_inf);
    gamma=1.4;
    T=288;
    R=287;
    C=(1-M.^2);
    Uinf=M*sqrt(gamma*T*R);

    %% Cell size
    dx=xdomain/cellx;
    dy=ydomain/celly;

    dxdy=dx/dy;
    dydx=dy/dx;
    %% Variables
    x=cellx;
    y=celly;
    total=x*y;

    %% Matrix A
    coef=sparse(total,total);
    tic
    k=1:1:y-2;
    for i=1:1:total
        flag=0;
        if ismember(i,1+k*x)==1                 %Middle-Left Cells
            coef(i,i)=(2*C*dydx-2*dxdy);
            coef(i,i+x)=dxdy;
            coef(i,i-x)=dxdy;
            flag=1;
        elseif ismember(i,2+k*x)                 %Middle-Left 2 column
            coef(i,i)=(C*dydx-2*dxdy);
            coef(i,i-1)=-3*C*dydx;
            coef(i,i+x)=dxdy;
            coef(i,i-x)=dxdy;
            flag=1;
        end
        if i==1                     %Bottom-Left Cell
            coef(i,i)=(2*C*dydx-dxdy);
            coef(i,i+x)=dxdy;
            flag=1;
        elseif i==2                 %Bottom-Left 2nd Cell
            coef(i,i)=C*dydx-dxdy;
            coef(i,i-1)=-3*C*dydx;
            coef(i,i+x)=dxdy;
            flag=1;
        elseif i==(x*(y-1)+1)       %Top-Left Cell modificado
            coef(i,i)=(2*C*dydx-dxdy);
            coef(i,i-x)=dxdy;
            flag=1;
        elseif i==(x*(y-1)+2)        %Top-Left 2nd Cell
            coef(i,i)=C*dydx-dxdy;
            coef(i,i-1)=-3*C*dydx;
            coef(i,i-x)=dxdy;
            flag=1;
        elseif (2<i && i<=x)         %Bottom Cells modificado, inclui bottom right
            coef(i,i)=(C*dydx-dxdy);
            coef(i,i-1)=-2*C*dydx;
            coef(i,i-2)=C*dydx;
            coef(i,i+x)=dxdy;
            flag=1;
        elseif (x*(y-1)+2<i && i<=x*y) %Top Cells and top right
            coef(i,i)=(C*dydx-dxdy);
            coef(i,i-1)=-2*C*dydx;
            coef(i,i-2)=C*dydx;
            coef(i,i-x)=dxdy;
            flag=1;
        end
        if flag==0                         %Middle Cells modificado
            coef(i,i)=(C*dydx-2*dxdy);
            coef(i,i-1)=-2*C*dydx;
            coef(i,i-2)=C*dydx;
            coef(i,i-x)=dxdy;
            coef(i,i+x)=dxdy;
        end
    end

    %% Matrix B
    value=sparse(total,1);

    for i=1:x
        if dx/2+(i-1)*dx>=airfoil_start && dx/2+(i-1)*dx<airfoil_start+c/2
            value(i,1)=slope*Uinf*dx;
        elseif dx/2+(i-1)*dx>airfoil_start+c/2 && dx/2+(i-1)*dx<=airfoil_start+c
            value(i,1)=-slope*Uinf*dx;
        end
    end

    %% Solution of Ax=B
    sol = coef\value;
    sol=full(sol);
    %% Potential Field Matrix
    phi=vec2mat(sol,x,y);
    %% u and v
    v=v_field(sol,x,y,dy,dx,Uinf,airfoil_start,c,slope);
    u=u_field(sol,x,y,dx);
    %% Calculations
    Utot =u+Uinf;                       %Total Velocity X
    Velocity=sqrt(Utot.^2+v.^2);        %Velocity Magnitude
    T0=T_inf*(1+(gamma-1)/2*M^2);       %Stagnation Temperature

    %Temperature Field
    T_field=(T0*R*gamma/(gamma-1)-0.5*Velocity.^2)*(gamma-1)/(gamma*R);

    %Mach Number Field
    Mach_local=Velocity./((gamma*R*T_field).^0.5);

    %Pressure Field
    P_field=p_inf*(T_field/T_inf).^(gamma/(gamma-1));

    %Stagnation Pressure Field
    P0_field=P_field.*(1+((gamma-1)/2)*Mach_local.^2).^(gamma/(gamma-1));

    %Density Field
    rho_field=rho_inf*(P_field/p_inf).^(1/gamma);

    %% Mach angle
    beta=shock_angle(phi,x,y)

    %% Plot of Potential Field
    X=linspace(dx,xdomain-dx,cellx);
    Y=linspace(dx,ydomain-dx,celly);
    surf(X,Y,phi,'EdgeColor','none');
    view(0,90)
    colorbar
    colormap('jet')
    %% Plot of Mach Flow Field
    figure(2)
    surf(X,Y,vec2mat(Mach_local,x),'EdgeColor','none');
    view(0,90);
    colorbar
    %% Plot of Density Field
    figure(3)
    surf(X,Y,vec2mat(rho_field,x),'EdgeColor','none');
    view(0,90);

    toc
    function [v] = v_field(sol,x,y,dy,dx,Uinf,airfoil_start,c,slope)
        total=x*y;
        v=zeros(total,1); %central differences
        for i=1:1:total
            if i<=x
                if dx/2+(i-1)*dx>=airfoil_start && dx/2+(i-1)*dx<airfoil_start+c/2 %Airfol ascent
                    v(i)=slope*Uinf;
                elseif dx/2+(i-1)*dx>airfoil_start+c/2 && dx/2+(i-1)*dx<=airfoil_start+c %Airfoil descent
                    v(i)=-slope*Uinf;
                else
                    v(i)=(sol(i+x)-sol(i))/dy;
                end
            elseif i>=1+x*(y-1)
                v(i)=(sol(i)-sol(i-x))/dy;

            elseif i>x && i<1+x*(y-1)
                v(i)=(sol(i+x)-sol(i-x))/(2*dy);
            end
        end
    end   
    function [u] = u_field(sol,x,y,dx)
        %% Compute U
        total=x*y;
        u=zeros(total,1); %central differences
        k=0:1:y-1;
        for i=1:1:total
            if ismember(i,1+k*x)==1
                u(i)= (sol(i))/dx;
            elseif ismember(i,x+k*x)==1
                u(i)= (sol(i)-sol(i-1))/dx;
            else                %SuperSonico
                u(i)= (sol(i)-sol(i-1))/dx;
            end
        end
    end
    function [beta]=shock_angle(phi,x,y)
        y_pos=[1:1:y/3]';
        for i=1:length(y_pos)
            for j=1:1:x
                if phi(y_pos(i),j)<-1
                    x_pos(i,1)=j;
                    break
                end
            end
        end
        m=(x_pos-x_pos(1))\y_pos;
        beta=atan(m)*360/(2*pi);
    end
end

function [sol,u,v,value,coef]=subsonic_field(M,xdomain,ydomain,cellx,celly,airfoil_start,c)   
%% Problem Constants
    p_inf=10^5;
    T_inf=288;
    R=287;
    rho_inf=p_inf/(R*T_inf);
    gamma=1.4;
    T=288;
    R=287;
    C=(1-M.^2);
    Uinf=M*sqrt(gamma*T*R);
    %% Cell size
    dx=xdomain/cellx;
    dy=ydomain/celly;

    dxdy=dx/dy;
    dydx=dy/dx;
    %% Variables
    x=cellx;
    y=celly;
    total=x*y;

    %% Matrix A
    coef=sparse(total,total);
    k=1:1:y-2;
    for i=1:1:total
        flag=0;
        if ismember(i,1+k*x)== 1                 %Middle-Left Cells
            coef(i,i)=-(3*C*dydx+2*dxdy);
            coef(i,i+1)=C*dydx;
            coef(i,i+x)=1*dxdy;
            coef(i,i-x)=1*dxdy;
            flag=1;
        elseif ismember(i,x+k*x)==1         %Middle-Right Cells
            coef(i,i)=-(3*C*dydx+2*dxdy);
            coef(i,i-1)=C*dydx;
            coef(i,i+x)=1*dxdy;
            coef(i,i-x)=1*dxdy;
            flag=1;
        end
        if i == total           %Top-Right Cell
            coef(i,i)=-(3*C*dydx+1*dxdy);
            coef(i,i-1)=C*dydx;
            coef(i,i-x)=1*dxdy;
            flag=1;
        elseif (i>x*(y-1)+1 && i<x*y)       %Top Cells
            coef(i,i)=-(2*C*dydx+1*dxdy);
            coef(i,i-1)=C*dydx;
            coef(i,i+1)=C*dydx;
            coef(i,i-x)=1*dxdy;
            flag=1;
        elseif i==(x*(y-1)+1)               %Top-Left Cell
            coef(i,i)=-(3*C*dydx+1*dxdy);
            coef(i,i+1)=C*dydx;
            coef(i,i-x)=1*dxdy;
            flag=1;
        elseif i==1             %Bottom-Left Cell
            coef(i,i)=-(3*C*dydx+1*dxdy);
            coef(i,i+1)=C*dydx;
            coef(i,i+x)=1*dxdy;
            flag=1;
        elseif i==x             %Bottom-Right
            coef(i,i)=-(3*C*dydx+1*dxdy);
            coef(i,i-1)=C*dydx;
            coef(i,i+x)=1*dxdy;
            flag=1;
        elseif (i>1 && i<x)             %Bottom Cells
            coef(i,i)=-(2*C*dydx+1*dxdy);
            coef(i,i-1)=C*dydx;
            coef(i,i+1)=C*dydx;
            coef(i,i+x)=1*dxdy;
            flag=1;
        end
        if flag==0                      %Middle Cells
            coef(i,i)=-(2*C*dydx+2*dxdy);
            coef(i,i-1)=C*dydx;
            coef(i,i+1)=C*dydx;
            coef(i,i-x)=1*dxdy;
            coef(i,i+x)=1*dxdy;
        end
    end

    %% Matrix B
    value=sparse(total,1);
    for i=1:x
        if dx/2+(i-1)*dx>=airfoil_start && dx/2+(i-1)*dx<airfoil_start+c/2
            value(i,1)=0.06*Uinf*dx;
        elseif dx/2+(i-1)*dx>= airfoil_start+c/2 && dx/2+(i-1)*dx<=airfoil_start+c
            value(i,1)=-0.06*Uinf*dx;
        end
    end
    %% Solution of Ax=B
    sol=coef\value;
    sol=full(sol);
    %% Potential Field Matrix
    phi=vec2mat(sol,x);
    %% Compute U e V
    u=u_field(sol,x,y,dx);
    v=v_field(sol,x,y,dy,dx,Uinf,airfoil_start,c);
    %% Calculations
    Utot =u+Uinf;                       %Total Velocity X
    Velocity=sqrt(Utot.^2+v.^2);        %Velocity Magnitude
    T0=T_inf*(1+(gamma-1)/2*M^2);       %Stagnation Temperature

    %Temperature Field
    T_field=(T0*R*gamma/(gamma-1)-0.5*Velocity.^2)*(gamma-1)/(gamma*R);

    %Mach Number Field
    Mach_local=Velocity./((gamma*R*T_field).^0.5);

    %Pressure Field
    P_field=p_inf*(T_field/T_inf).^(gamma/(gamma-1));

    %Stagnation Pressure Field
    P0_field=P_field.*(1+((gamma-1)/2)*Mach_local.^2).^(gamma/(gamma-1));

    %Density Field
    rho_field=rho_inf*(P_field/p_inf).^(1/gamma);



    %% Plot of Potential Field
    figure(1)
    X=linspace(dx,xdomain-dx,cellx);
    Y=linspace(dx,ydomain-dx,celly);
    surf(X,Y,phi,'EdgeColor','none','FaceColor','interp');
    view(0,90)
    colorbar
    colormap('hot')
    %% Plot of Mach Flow Field
    figure(2)
    surf(X,Y,vec2mat(Mach_local,x),'EdgeColor','none');
    view(0,90);
    colorbar
    %% Plot of Density Field
    figure(3)
    surf(X,Y,vec2mat(rho_field,x),'EdgeColor','none');
    view(0,90);


    function [u] = u_field(sol,x,y,dx,K)
        %% Compute U
        total=x*y;
        u=zeros(total,1); %central differences
        k=0:1:y-1;
        for i=1:1:total
            if ismember(i,1+k*x)==1
                u(i)= (sol(i+1)-sol(i))/dx;
            elseif ismember(i,x+k*x)==1
                u(i)= (sol(i)-sol(i-1))/dx;
            else
                u(i)= (sol(i+1)-sol(i-1))/(2*dx);
            end
        end
    end


    function [v] = v_field(sol,x,y,dy,dx,Uinf,airfoil_start,c)
        total=x*y;
        v=zeros(total,1); %central differences
        for i=1:1:total
            if i<=x
                if dx/2+(i-1)*dx>=airfoil_start && dx/2+(i-1)*dx<airfoil_start+c/2 %Airfol ascent
                    v(i)=0.06*Uinf;
                elseif dx/2+(i-1)*dx>airfoil_start+c/2 && dx/2+(i-1)*dx<=airfoil_start+c %Airfoil descent
                    v(i)=-0.06*Uinf;
                else
                    v(i)=(sol(i+x)-sol(i))/dy;
                end
            elseif i>=1+x*(y-1)
                v(i)=(sol(i)-sol(i-x))/dy;

            elseif i>x && i<1+x*(y-1)
                v(i)=(sol(i+x)-sol(i-x))/(2*dy);
            end
        end
    end
end
    
%% Jacobi's Iterative Method
function x0=Jacobi_Method(A,b)
%This function uses Jacobi's iterative method to solve the variable in a
%system of equations

%% -----INPUT-----

%Input A-matrix (square)
% 
% A = [ ];
% 
% %Input corresponding b-matrix (column)
% 
% b = [];

%Specify number of iterations to run

itr = 3;

%Unless otherwise specified, the intitial guess x_0 will be all zeros

x0 = zeros(1,length(A));

%% -----SOLUTION-----

n = length(A);

%To keep track of the results, we define a matrix (x) as

x = [x0;zeros(itr,n)];
    
%The iterative method below will run for the specified number of iteration
%above, calculate, and tabulate the values of x^k.

for k = 1:itr
    for i = 1:n
        
        %We define a variable "sigma" that will be used to sum the values
        %that are known (not on the main diagonal of the matrix) for each
        %row of the A-matrix (performing the calculations for each
        %equation)
        sigma = 0;
        
        for j = 1:n
            
            %Because the coefficient along the main diagonal is paired with
            %the x_i we are solving for, it is omitted from the sum
            if i~=j
                sigma = sigma + A(i,j)*x(k,j);
            end
            
        end
        
        %Lastly, the x_i is calculated from the equation and recorded under
        %its respective iteration
        x(k+1,i) = (b(i)-sigma)/A(i,i);
        
    end
    
end

%% -----OUTPUT-----
%Output created for 3x3 A-matrix

%For output purposes, 'k' will be defined to show the iteration number
k = [0:itr]';
end
