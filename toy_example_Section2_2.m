clear; clc;  
% 设置全局参数数组
global Kt  Nt Nx dt T_end flag_eigPA theta alpha flag_GMRES
alpha=5e-3;
theta=1/2;
flag_eigPA=1;
flag_GMRES=1;
Nx=3;
% 时间离散参数
dt = 0.025;          % 时间步长
T_end =5;    

Nt = ceil(T_end / dt);  % 时间步数
time = 0:dt:T_end;  % 时间网格
It=speye(Nt);
 
% 生成三角形网格
 
 Kt=zeros(Nx,Nx,Nt);
 for n=1:Nt
     Kt(:,:,n)= Kt_fun(time(n)+dt*theta);
 end 
Err0=zeros(1,Nt);
Err2=zeros(1,Nt);
Err3=zeros(1,Nt);
[barK0,b0]=getKb(0);
[barK2,b2]=getKb(2);
[barK3,b3]=getKb(3);
for n=1:Nt
    Err0(n)=norm(Kt(:,:,n)-barK0,'fro');
    Err2(n)=norm(Kt(:,:,n)-b2(n)*barK2,'fro');
    Err3(n)=norm(Kt(:,:,n)-b3(n)*barK3,'fro');
end
 
figure(1);
semilogy(1:Nt,Err0,'k-.',1:Nt,Err2,'r-.',1:Nt,Err3,'b-.','LineWidth',1,'markersize',12);
%plot(1:Nt,Err2,'r-.',1:Nt,Err3,'b-.','LineWidth',1,'markersize',12);

set(gca,'FontSize',15);
xlabel('$n$ (index of time point)','FontSize',21,'Interpreter','latex');
ylabel('$\|K_n-d_n\overline{K}\|_{\rm F}$','FontSize',21,'Interpreter','latex');
leg=legend('Method-1', 'Method-2','Method-3');
set(leg,'fontsize',16);
xlim([1,Nt]);
% if Nt==750
%     set(gca,'xtick',[1,150:150:Nt]);
% elseif Nt==500
%     set(gca,'xtick',[1,100:100:Nt]);
% else
%     set(gca,'xtick',[1,50:50:Nt]);
% end
%title(['$N_t=',num2str(Nt),'$'],'Interpreter','latex','FontSize',21);
%ylim([0,8]);
grid on;

if flag_eigPA==1
    e=ones(Nt,1);
    B1=spdiags([-e e], [-1,0], Nt, Nt);
    B2=spdiags([(1-theta)*e theta*e], [-1,0], Nt, Nt);
    C1=B1;
    C1(1,Nt)=alpha*C1(2,1);
    C2=B2;
    C2(1,Nt)=alpha*C2(2,1); 
    It=speye(Nt);Ix=speye(Nx);
    calK=sparse(Nt*Nx,Nt*Nx);
    for n=1:Nt
         calK((n-1)*Nx+1:n*Nx,(n-1)*Nx+1:n*Nx)=Kt(:,:,n);
    end
    calA=kron(B1,Ix)+dt*calK*kron(B2,Ix);
    [barK0,b0]=getKb(0);
     calP0=kron(C1,Ix)+dt*kron(diag(b0)*C2,barK0);
    [barK2,b2]=getKb(2);
     calP2=kron(C1,Ix)+dt*kron(diag(b2)*C2,barK2);
    [barK3,b3]=getKb(3);
     calP3=kron(C1,Ix)+dt*kron(diag(b3)*C2,barK3);
    fprintf('计算预处理矩阵特征值 for Method-1...\n');
    eig_invP0A=eig(full(calP0\calA));
    fprintf('计算预处理矩阵特征值 for Method-2...\n');
    eig_invP2A=eig(full(calP2\calA));
    fprintf('计算预处理矩阵特征值 for Method-3...\n');
    eig_invP3A=eig(full(calP3\calA));
    figure(2);
    plot(real(eig_invP0A),imag(eig_invP0A),'ko');
    set(gca,'FontSize',15);
    xlabel('Re$(\mathcal{P}_{\alpha}^{-1}\mathcal{A})$','FontSize',21,'Interpreter','latex');
    ylabel('Im$(\mathcal{P}_{\alpha}^{-1}\mathcal{A})$','FontSize',21,'Interpreter','latex');
    title('Method-1','FontSize',21);
    xlim([min(real(eig_invP0A))-0.1/2,max(real(eig_invP0A))+0.1/2]);
    ylim([min(imag(eig_invP0A))-0.1/2,max(imag(eig_invP0A))+0.1/2]);

    figure(3);
    plot(real(eig_invP2A),imag(eig_invP2A),'r*');
    set(gca,'FontSize',15);
    xlabel('Re$(\mathcal{P}_{\alpha}^{-1}\mathcal{A})$','FontSize',21,'Interpreter','latex');
    ylabel('Im$(\mathcal{P}_{\alpha}^{-1}\mathcal{A})$','FontSize',21,'Interpreter','latex');
    title('Method-2','FontSize',21);
    xlim([min(real(eig_invP0A))-0.1/2,max(real(eig_invP0A))+0.1/2]);
    ylim([min(imag(eig_invP0A))-0.1/2,max(imag(eig_invP0A))+0.1/2]);
    figure(4);
    plot(real(eig_invP3A),imag(eig_invP3A),'bs');
    set(gca,'FontSize',15);
    xlabel('Re$(\mathcal{P}_{\alpha}^{-1}\mathcal{A})$','FontSize',21,'Interpreter','latex');
    ylabel('Im$(\mathcal{P}_{\alpha}^{-1}\mathcal{A})$','FontSize',21,'Interpreter','latex');
    title('Method-3','FontSize',21);
    xlim([min(real(eig_invP0A))-0.1/2,max(real(eig_invP0A))+0.1/2]);
    ylim([min(imag(eig_invP0A))-0.1/2,max(imag(eig_invP0A))+0.1/2]);
    if flag_GMRES==1
        b=zeros(Nx*Nt,1);
        U0=[0;0;0];
        source_term=zeros(Nx,Nt);
        for n=1:Nt
            tn=(n+theta)*dt;
            source_term(:,n)=[1-tn;5*tn^0.5;-0.5*tn];
            %source_term(:,n)=0*random('Uniform',-1,1,Nx,1);
        end
        for n=1:Nt
            if n==1
                b((n-1)*Nx+1:n*Nx)=(Ix - dt*(1-theta) * Kt(:,:,n)) * U0 +dt* source_term(:,n);
            else
                b((n-1)*Nx+1:n*Nx)=dt*source_term(:,n);
            end
        end
        U_Ini=kron(ones(Nt,1),1+U0);
        %b=random('Uniform',-1,1,Nx*Nt,1);
        [U_ref0,FLAG,RELRES,ITER,RESVEC0]=gmres(calA,b,20,1e-10,50,...
                       calP0,[],U_Ini);
        [U_ref2,FLAG,RELRES,ITER,RESVEC2]=gmres(calA,b,20,1e-10,50,...
                       calP2,[],U_Ini);
        [U_ref3,FLAG,RELRES,ITER,RESVEC3]=gmres(calA,b,20,1e-10,50,...
                           calP3,[],U_Ini);
        figure(5);
        semilogy(0:length(RESVEC0)-1,RESVEC0/max(RESVEC0),'k-.o',...
            0:length(RESVEC2)-1,RESVEC2/max(RESVEC2),'r-.*',...
            0:length(RESVEC3)-1,RESVEC3/max(RESVEC3),'b-.s','MarkerSize',12,'LineWidth',1);
        shg;
        set(gca, 'FontSize',15);
        leg=legend('Method-1', 'Method-2','Method-3');
        set(leg,'fontsize',16);
        xlabel('GMRES iteration number','FontSize',21);
        ylabel('Residual','FontSize',21);
        xlim([0,11]);
        set(gca,'xtick',0:11);
        set(gca,'ytick',10.^(-10:0));
        ylim([1e-10,1.5]);
        title(['$N_t=',num2str(Nt),'$'],'Interpreter','latex','FontSize',21);
    end
end
function [barK,b]=getKb(solver_P)
global Nt Nx Kt
if solver_P==0
    barK=zeros(Nx,Nx);
    for n=1:Nt
        barK=barK+Kt(:,:,n)/Nt;
    end
    %val=(kron(B1,Ix)+dt*kron(B2,barK))\du;
    b=ones(Nt,1);
elseif solver_P==2
     barK=zeros(Nx,Nx);
    for n=1:Nt
        barK=barK+Kt(:,:,n)/Nt;
    end
    denominator=trace(transpose(barK)*barK);
    for n=1:Nt
        b(n)=trace(transpose(Kt(:,:,n))*barK)/denominator;
    end
else
    A1 = Kt(:,:,2);
    [m, ~] = size(A1);
    vecs = zeros(m^2, Nt);
    for n = 1:Nt
        An = Kt(:,:,n);   % 获取第n个块矩阵
        vecs(:, n) = An(:);  % 向量化
    end
    % 计算Gram矩阵
    G = vecs' * vecs;
    [V, eigvals] = eigs(G, 1,'la');
    u = (V(:, 1));  lambda = eigvals(1);  
    u = u / norm(u,2);
    barK = sparse(m, m);
    for n = 1:Nt
        barK = barK + u(n) * Kt(:,:,n);
    end
    b = u;
    % if sum(b < 0) > Nt/2
    %     b = -b;
    %     barK = -barK;
    % end
end
end


function val=Kt_fun(t)
global T_end
V=[1,-1,0;
     t,t,0; 
     t,t,1];
%D=diag([0;2/2-sin(1*sqrt(t))^2;0+t^2])*2;
D=diag([0;1-sin(0.5*t);t^2])*2;
val=V*D*inv(V); 
end