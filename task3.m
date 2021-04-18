%driven lid cavity flow
%% setups

clc;
clear;
xmax=1;
ymax=1;
Re=0.1;%Reynolds number
Nx=41;%Number of nodes in the horizontal direction
Ny=41;%Number of nodes in the vertical direction

U_inf=1;%Lid velocity

DeltaX=xmax/(Nx-1);
DeltaY=ymax/(Ny-1);
if Re == 0.1
    DeltaT=0.00001;
elseif Re == 1
    DeltaT=0.0001;
end
x=linspace(0,xmax,Nx);
y=linspace(0,ymax,Ny);
[X,Y]=meshgrid(x,y);

nt=0.1;%Total simulation time
u(:,:,1)=zeros(Ny,Nx);
v(:,:,1)=zeros(Ny,Nx);
for i=2:(Nx-1)
    u(Ny,i,1)=U_inf;%Initialization of moving lid
end

%Initialization of the velocity field
p(:,:,1)=zeros(Ny,Nx);%Initialization of the pressure field
Converginglevel=0.001;%Converging level for pressure-Poisson equations
Converginglevel1=0.000001;%Converging level for x-velocity
Converginglevel2=0.000001;%Converging level for y-velocity
uint=zeros(Ny,Nx);
vint=zeros(Ny,Nx);
CE=zeros(Ny,Nx);%Matrix to store continuity errors.

%% Iteration
tic;
for n=1:round(nt/DeltaT)
    %Predict Velocity
    for i=2:(Nx-1)
        for j=2:(Ny-1)
            Cx=u(j,i,n)/(2*DeltaX)*(u(j,i+1,n)-u(j,i-1,n))+v(j,i,n)/(2*DeltaY)*(u(j+1,i,n)-u(j-1,i,n));%Horizontal convective term
            Cy=u(j,i,n)/(2*DeltaX)*(v(j,i+1,n)-v(j,i-1,n))+v(j,i,n)/(2*DeltaY)*(v(j+1,i,n)-v(j-1,i,n));%Vertical convective term
            Dx=(u(j,i+1,n)-2*u(j,i,n)+u(j,i-1,n))/(Re*DeltaX^2)+(u(j+1,i,n)-2*u(j,i,n)+u(j-1,i,n))/(Re*DeltaY^2);%Horizontal diffusive term
            Dy=(v(j,i+1,n)-2*v(j,i,n)+v(j,i-1,n))/(Re*DeltaX^2)+(v(j+1,i,n)-2*v(j,i,n)+v(j-1,i,n))/(Re*DeltaY^2);%Vetical diffusive term
            uint(j,i)=u(j,i,n)+DeltaT*(-1*Cx+Dx);
            vint(j,i)=v(j,i,n)+DeltaT*(-1*Cy+Dy);
        end
    end
    

    for i=2:(Nx-1)
        for j=2:(Ny-1)
            CE(j,i)=(uint(j,i+1)-uint(j,i-1))/(2*DeltaX*DeltaT)+(vint(j+1,i)-vint(j-1,i))/(2*DeltaY*DeltaT);
            %Continuty error
        end
    end
    
    
    %Using SOR to solve the Pressure-Poisson Equations.
    omega=1.4;
    p(:,:,n+1)=p(:,:,n);
    
    
    for iii=1:200
        for i=2:(Nx-1)
            p(1,i,n+1)=p(2,i,n+1);
            p(Ny,i,n+1)=p(Ny-1,i,n+1);
        end
        for j=2:(Ny-1)
            p(j,1,n+1)=p(j,2,n+1);
            p(j,Nx,n+1)=p(j,Nx-1,n+1);
        end
        for i=2:(Nx-1)
            for j=2:(Ny-1)
                R=((p(j,i+1,n+1)+p(j,i-1,n+1))/(DeltaX^2)+(p(j+1,i,n+1)+p(j-1,i,n+1))/(DeltaY^2)-CE(j,i))/(2/(DeltaX^2)+2/(DeltaY^2));
                p(j,i,n+1)=omega*R+(1-omega)*p(j,i,n+1);
            end
        end
        
        residual=0;%Calculation error residual for pressure-Poisson equations
        for i=2:(Nx-1)
            for j=2:(Ny-1)
                e=CE(j,i)-(p(j,i+1,n+1)-2*p(j,i,n+1)+p(j,i-1,n+1))/(DeltaX^2)-(p(j+1,i,n+1)-2*p(j,i,n+1)+p(j-1,i,n+1))/(DeltaY^2);
                residual=residual+e^2;
            end
        end
        residual=sqrt(residual/((Nx-1)*(Ny-1)));
        if residual<Converginglevel
            break
        end
    end
    
    for i=2:(Nx-1)
        p(1,i,n+1)=p(2,i,n+1);
        p(Ny,i,n+1)=p(Ny-1,i,n+1);
    end
    for j=2:(Ny-1)
        p(j,1,n+1)=p(j,2,n+1);
        p(j,Nx,n+1)=p(j,Nx-1,n+1);
    end
    
    residual1=0;%Initialization of residual in x-velocity
    residual2=0;%Initialization of residual in y-velocity
    
    %Correction of grid node velocities
    for i=2:(Nx-1)
        for j=2:(Ny-1)
            u(j,i,n+1)=uint(j,i)-DeltaT*(p(j,i+1,n+1)-p(j,i-1,n+1))/(2*DeltaX);
            v(j,i,n+1)=vint(j,i)-DeltaT*(p(j+1,i,n+1)-p(j-1,i,n+1))/(2*DeltaY);
        end
    end    
   
    %Boundary Condition
    u(Ny,2:Nx-1,n+1)=U_inf;
    u(1,2:Nx-1,n+1)=0;
    v(Ny,2:Nx-1,n+1)=0;
    v(1,2:Nx-1,n+1)=0;
    u(2:Ny-1,1,n+1)=0;
    u(2:Ny-1,Nx,n+1)=0;
    v(2:Ny-1,1,n+1)=0;
    v(2:Ny-1,Nx,n+1)=0;

end
time_elapsed = toc;
%% Post Process
time_stamps = [0.02, 0.04, 0.06, 0.08, 0.1];
time_step_values = time_stamps*n/nt + 1;

%Streamlines

startx1=[];
startx2=[];
startx3=[];
starty1=[];
starty2=[];
starty3=[];
FigHandle = figure('Position', [100, 40, 800, 600]);
Factor=1;%Coefficient of pathline skipping
startx1=0:Factor*DeltaX:xmax-2*Factor*DeltaX;
startx2=Factor*DeltaX:Factor*DeltaX:xmax-Factor*DeltaX;
startx3=DeltaX*Factor*2:Factor*DeltaX:xmax;
Factor=2;
for i=0:Factor*DeltaY:ymax
    starty1(round(1+i/DeltaY/Factor),:)=i*ones(size(startx1));
end
for i=0:2*Factor*DeltaY:ymax
    starty2(round(1+i/DeltaY/Factor),:)=i*ones(size(startx2));
end
for i=0:3*Factor*DeltaY:ymax
    starty3(round(1+i/DeltaY/Factor),:)=i*ones(size(startx3));
end

%quiver([DeltaX/2:DeltaX:xmax-DeltaX/2],[DeltaY/2:DeltaY:ymax-DeltaY/2],u(2:Ny-1,2:Nx-1,n),v(2:Ny-1,2:Nx-1,n));
figure(1)
tiledlayout(2,3);
nexttile
for i=0:Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(1)),startx1,starty1(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','red');
end
for i=0:2*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(1)),startx2,starty2(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','blue');
end
for i=0:3*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(1)),startx3,starty3(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','green');
end
axis([0 xmax 0 ymax]);
axis square;
title('At t = 0.02')
xlabel('x(m)');
ylabel('y(m)');

nexttile
for i=0:Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(2)),startx1,starty1(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','red');
end
for i=0:2*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(2)),startx2,starty2(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','blue');
end
for i=0:3*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(2)),startx3,starty3(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','green');
end
axis([0 xmax 0 ymax]);
axis square;
title('At t = 0.04')
xlabel('x(m)');
ylabel('y(m)');

nexttile
for i=0:Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(3)),startx1,starty1(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','red');
end
for i=0:2*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(3)),startx2,starty2(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','blue');
end
for i=0:3*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(3)),startx3,starty3(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','green');
end
axis([0 xmax 0 ymax]);
axis square;
title('At t = 0.06')
xlabel('x(m)');
ylabel('y(m)');

nexttile
for i=0:Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(4)),startx1,starty1(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','red');
end
for i=0:2*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(4)),startx2,starty2(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','blue');
end
for i=0:3*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(4)),startx3,starty3(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','green');
end
axis([0 xmax 0 ymax]);
axis square;
title('At t = 0.08')
xlabel('x(m)');
ylabel('y(m)');

nexttile
for i=0:Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(5)),startx1,starty1(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','red');
end
for i=0:2*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(5)),startx2,starty2(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','blue');
end
for i=0:3*Factor*DeltaY:ymax
    h=streamline(X,Y,u(:,:,n),v(:,:,time_step_values(5)),startx3,starty3(round(1+i/DeltaY/Factor),:),[0.1,500]);set(h,'Color','green');
end
axis([0 xmax 0 ymax]);
axis square;
title('At t = 0.1')
xlabel('x(m)');
ylabel('y(m)');

%Contour of Pressure
figure(3)
tiledlayout(2,3);
nexttile
contourf(X,Y,p(:,:,time_step_values(1)),20,'ShowText','on')
title('Pressure at t = 0.02')
xlabel('x(m)');
ylabel('y(m)');
nexttile
contourf(X,Y,p(:,:,time_step_values(2)),20,'ShowText','on')
title('Pressure at t = 0.04')
xlabel('x(m)');
ylabel('y(m)');
nexttile
contourf(X,Y,p(:,:,time_step_values(3)),20,'ShowText','on')
title('Pressure at t = 0.06')
xlabel('x(m)');
ylabel('y(m)');
nexttile
contourf(X,Y,p(:,:,time_step_values(4)),20,'ShowText','on')
title('Pressure at t = 0.08')
xlabel('x(m)');
ylabel('y(m)');
nexttile
contourf(X,Y,p(:,:,time_step_values(5)),20,'ShowText','on')
title('Pressure at t = 0.1')
xlabel('x(m)');
ylabel('y(m)');

h=colorbar;
title(h,'pressure');
h.Layout.Tile = 'east';

%u-velocity profile on vertical center line
figure(4)
plot(0:DeltaY:ymax,u(:,round((1+Nx)/2),n));
title('u-velocity profile on vertical center line')
xlabel('y(m)');
ylabel('u at x=0.5');

%v-velocity profile on horizontal center line
figure(5)
plot(0:DeltaX:xmax,v(:,round((1+Nx)/2),n));
title('v-velocity profile on vertical center line')
xlabel('y(m)');
ylabel('v at x=0.5');