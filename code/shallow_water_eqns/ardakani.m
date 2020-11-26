% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Shallow-water sloshing in vessels undergoing prescribed rigid-body motion
% in two dimensions.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Written by H. Alemi Ardakani & T.J. Bridges,
% Department of Mathematics, University of Surrey, UK.


clear
clc
close all;

tic

format long

g=9.81;
L=0.50;
N=100;
dx=L/N;

dt=1E-2;
d1=-L/2;
d2=0.0;

h0=0.08;
eps_p=3.0*pi/180;                 % pitching force amplitude
eps_s=0.02;                       % surging force amplitude
eps_h=0.015;                      % heaving force amplitude
omega_1=(pi/L)*sqrt(g*h0);
omega_2=(2*pi/L)*sqrt(g*h0);
omega_n=0.9*omega_1;              % pitch motion frequency
omega_s=0.8*omega_1;              % surge motion frequency
omega_h=0.8*omega_1;              % heave motion frequency


for i=1:N+1
    X(i)=(i-1)*dx;
end

EH(1:N+1)=h0;
Ufirst(1:N+1)=0.0;

t(1)=0;
T=t(1);

% Plotting the sloshing plate at t=0.0 sec:
theta=eps_p*sin(omega_n*T);
q_1=eps_s*sin(omega_s*T);
q_2=eps_h*sin(omega_h*T);

Q=[cos(theta) -sin(theta);sin(theta) cos(theta)];

v_1=Q*[d1;d2]+[q_1;q_2];
v_2=Q*[d1+L;d2]+[q_1;q_2];
v_3=Q*[d1;d2+(L/2.0)]+[q_1;q_2];
v_4=Q*[d1+L;d2+(L/2.0)]+[q_1;q_2];

box_1=[-L;-L/2.0];
box_2=[L;-L/2.0];
box_3=[-L;L];
box_4=[L;L];

v_r=[q_1;q_2];

for jj=1:N+1
    H_p{1,jj}=(Q*[X(jj)+d1;EH(jj)+d2])+[q_1;q_2];
end
Hp=cell2mat(H_p);
figure()

set(gcf,'Color',[1,1,1])
line([box_1(1) box_2(1)],[box_1(2) box_2(2)],'Color','k','LineWidth',2.5);
line([box_1(1) box_3(1)],[box_1(2) box_3(2)],'Color','k','LineWidth',2.5);
line([box_3(1) box_4(1)],[box_3(2) box_4(2)],'Color','k','LineWidth',2.5);
line([box_2(1) box_4(1)],[box_2(2) box_4(2)],'Color','k','LineWidth',2.5);

line([v_1(1) v_2(1)],[v_1(2) v_2(2)],'Color','k','LineWidth',2.5);
line([v_1(1) v_3(1)],[v_1(2) v_3(2)],'Color','k','LineWidth',2.5);
line([v_3(1) v_4(1)],[v_3(2) v_4(2)],'Color','k','LineWidth',2.5);
line([v_2(1) v_4(1)],[v_2(2) v_4(2)],'Color','k','LineWidth',2.5);
line([box_1(1) box_2(1)],[0 0],'Color','r','LineWidth',2.5);

hold on
plot(v_r(1),v_r(2),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
plot(Hp(1,:),Hp(2,:),'b','LineWidth',2.5)
hold off
axis off
set(gca,'DataAspectRatio',[1 1 1])
Image(1)=getframe;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

for j=1:1000
    
    t(j+1)=j*dt;
    T=t(j+1);
    
    theta=eps_p*sin(omega_n*T); % counterclockwise pitch rotation
    bigomega=eps_p*omega_n*cos(omega_n*T);
    bigomegadot=-eps_p*(omega_n^2)*sin(omega_n*T);
    
    q_1=eps_s*sin(omega_s*T);
    qddot_1=-eps_s*(omega_s^2)*sin(omega_s*T); % surge forcing
    
    q_2=eps_h*sin(omega_h*T);
    qddot_2=-eps_h*(omega_h^2)*sin(omega_h*T); % heave forcing
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
    alpha=(g*cos(theta))+(bigomegadot.*(X+d1))-((bigomega^2)*d2)...
        -(qddot_1*sin(theta))+(qddot_2*cos(theta));
    
    beta=(-g*sin(theta))+(bigomegadot*d2)+((bigomega^2).*(X+d1))...
        -(qddot_1*cos(theta))-(qddot_2*sin(theta));
    
    A=cell(N+1);
    b=cell(N+1,1);
    Uzero=Ufirst;
    Hzero=EH;
    
    
    id=1;
    qpp=1;
    while qpp>1E-8
        
        alphatilde=alpha-((bigomega^2).*Hzero)+(2*bigomega.*Uzero);
        
        
        A{1,1}=[1 0;0 1];
        A{1,2}=(Hzero(1)*dt/dx)*[0 1;0 0];
        b{1,1}=[1 0;0 0]*[EH(1);Ufirst(1)];
        
        for i=3:N+1
            A{1,i}=zeros(2);
        end
        
        for ii=2:N
            
            if ii>2
                for i=1:ii-2
                    A{ii,i}=zeros(2);
                end
            end
            
            Atilde=(dt/(2.0*dx))*[alphatilde(ii) (Uzero(ii)...
                +(2.0*bigomega*Hzero(ii)));Uzero(ii) Hzero(ii)];
            
            A{ii,ii-1}=-Atilde;
            A{ii,ii}=[-bigomegadot*dt 1;1 0];
            A{ii,ii+1}=Atilde;
            
            for i=ii+2:N+1
                A{ii,i}=zeros(2);
            end
            b{ii,1}=([0 1;1 0]*[EH(ii);Ufirst(ii)])+(dt*beta(ii)*[1;0]);
        end
        
        for i=1:N-1
            A{N+1,i}=zeros(2);
        end
        
        A{N+1,N}=(-Hzero(N+1)*dt/dx)*[0 1;0 0];
        A{N+1,N+1}=[1 0;0 1];
        b{N+1,1}=[1 0;0 0]*[EH(N+1);Ufirst(N+1)];
        
        
        AA=cell2mat(A);
        bb=cell2mat(b);
        
        x=AA\bb;
        
        hsec=x(1:2:end);
        usec=x(2:2:end);
        Usec=usec';
        Hsec=hsec';
        
        qpp=max(abs(Hsec-Hzero)+abs(Usec-Uzero));
        
        id=id+1;
        
        Uzero=Usec;
        Hzero=Hsec;
        idd(j)=id;
        
    end % while
    
    Ufirst=Usec;
    EH=Hsec;
    
    % Plotting the sloshing plate at t(j+1)=j*dt:
    Q=[cos(theta) -sin(theta);sin(theta) cos(theta)];
    
    v_1=Q*[d1;d2]+[q_1;q_2];
    v_2=Q*[d1+L;d2]+[q_1;q_2];
    v_3=Q*[d1;d2+(L/2.0)]+[q_1;q_2];
    v_4=Q*[d1+L;d2+(L/2.0)]+[q_1;q_2];
    
    v_r=[q_1;q_2];
    
    l_1x=linspace(v_1(1),v_2(1),101);
    l_1y=linspace(v_1(2),v_2(2),101);
    l_2x=linspace(v_1(1),v_3(1),101);
    l_2y=linspace(v_1(2),v_3(2),101);
    l_3x=linspace(v_2(1),v_4(1),101);
    l_3y=linspace(v_2(2),v_4(2),101);
    l_4x=linspace(v_3(1),v_4(1),101);
    l_4y=linspace(v_3(2),v_4(2),101);
    
    for jj=1:N+1
        H_p{1,jj}=(Q*[X(jj)+d1;EH(jj)+d2])+[q_1;q_2];
    end
    Hp=cell2mat(H_p);
    
    set(gcf,'Color',[1,1,1])
    plot(l_1x,l_1y,'Color','k','LineWidth',2.5);
    hold on
    plot(l_2x,l_2y,'Color','k','LineWidth',2.5);
    plot(l_3x,l_3y,'Color','k','LineWidth',2.5);
    plot(l_4x,l_4y,'Color','k','LineWidth',2.5);
    plot(v_r(1),v_r(2),'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',8);
    plot(Hp(1,:),Hp(2,:),'b','LineWidth',2.5)
    hold off
    
    line([box_1(1) box_2(1)],[box_1(2) box_2(2)],'Color','k','LineWidth',2.5);
    line([box_1(1) box_3(1)],[box_1(2) box_3(2)],'Color','k','LineWidth',2.5);
    line([box_3(1) box_4(1)],[box_3(2) box_4(2)],'Color','k','LineWidth',2.5);
    line([box_2(1) box_4(1)],[box_2(2) box_4(2)],'Color','k','LineWidth',2.5);
    line([box_1(1) box_2(1)],[0 0],'Color','r','LineWidth',2.5);
    
    axis off
    set(gca,'DataAspectRatio',[1 1 1])
    Image(j+1)=getframe;
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
    % % % % % % % % % % Cr & Fr numbers:
    Cr=(dt/dx).*abs(Usec);
    Frs=(Usec.^2)/(g.*Hsec);
    Frsright=cos(theta)+((bigomegadot/g).*(X+d1))-(((bigomega^2)/g).*(Hsec+d2))...
        -((qddot_1/g)*sin(theta))+((qddot_2/g)*cos(theta));
    
    Crmax(j+1)=max(Cr);
    Frmax(j+1)=max(Frs);
    Frmaxright(j+1)=max(Frsright);
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    
end % j

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% movie(Image,1,50)
% map=colormap
% mpgwrite(Image,map,'Sloshing.mpg')
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

toc
elapsed_time=toc;