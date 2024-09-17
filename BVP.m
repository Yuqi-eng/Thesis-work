clc
close all

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 32)
set(0,'defaultaxeslinewidth',2)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',3)
set(0,'defaultTextInterpreter','latex')
%% Numerical set up
par.K = 30;                    % Number of spatial grid cells
par.L = 1;                      % Domain length
x = linspace(0,par.L,par.K);    % Discretise spatial domain

%% Plot figure
%%% Create initial guess
% solinit = bvpinit(x, @lowguess);
solinit = bvpinit(x, @highguess);

%%% Solve
% options = bvpset(Stats="on");
% sol = bvp5c(@bvpfcn, @bcfcn, solinit, options);
% 
% %%% reshape output
% n = sol.y(1,:);
% p = sol.y(3,:);
% u = sol.y(5,:);
% sig = sol.y(6,:);
% % plot(sol.x,[n;p;u]);
% figure('Units','normalized','Position',[0 0 0.3 0.3])
% subplot(1,4,1)
% plot(sol.x,n)
% % ylim([0,4])
% title('$n(t,x)$')
% axis square
% subplot(1,4,2)
% plot(sol.x,p)
% title('$\rho(t,x)$')
% axis square
% subplot(1,4,3)
% plot(sol.x,u)
% title('$u(t,x)$')
% axis square
% subplot(1,4,4)
% plot(sol.x,sig)
% title('$\sigma(t,x)$')
% axis square
% pic_name = ['BVP_res' '.png'];
% saveas(gcf,pic_name);

% Generate results for comparison with analytic solution
solinit1 = bvpinit(x, @lowguess);
solinit2 = bvpinit(x, @highguess);
options = bvpset(Stats="on");
sol1 = bvp5c(@bvpfcn, @bcfcn, solinit1, options);
sol2 = bvp5c(@bvpfcn, @bcfcn, solinit2, options);
[~,fullsol1] = Simulations_1DMix(2,2,0,par.K,4,10);
[x,fullsol2] = Simulations_1DMix(2,2,4,par.K,4,10);

%%% reshape output
n1 = sol1.y(1,:);
f1 = ones(length(x));
fulln1 = fullsol1(1:par.K,end);
n2 = sol2.y(1,:);
f2 = 5-4*(exp(x)+exp(-x))/(exp(1)+exp(-1));
fulln2 = fullsol2(1:par.K,end);
figure('Units','normalized','Position',[0 0 0.35 0.5])
hold on
pf1 = plot(x,f1,"blue");
pf2 = plot(x,f2,"blue");
pa1 = plot(sol1.x, n1,"+","LineWidth",2,"Color","red","MarkerSize",12);
pa2 = plot(sol2.x, n2,"+","LineWidth",2,"Color","red","MarkerSize",12);
pb1 = plot(x,fulln1,"o","LineWidth",2,"Color","green","MarkerSize",12);
pb2 = plot(x,fulln2,"o","LineWidth",2,"Color","green","MarkerSize",12);
ylim([0,4])
xlabel('Distance, $x$') 
ylabel('Cell density, $n(x)$')
hold off
legend([pf1(1), pf2(1), pa1(1), pa2(1) pb1(1) pb2(1)], {'analytic solution: normal', 'analytic solution: elevated', 'BVP: normal', 'BVP: elevated', 'PDE: normal', 'PDE: elevated'},'FontSize',18);
pic_name = ['BVP_compare' '.png'];
saveas(gcf,pic_name);
%% Main function implementing the model
function f = bvpfcn(~,y)
    %%% Parameter values
    a = 4;
    D2 = 0.1;
    d2 = 1;
    E1 = 100;
    taup2 = 10;
    
    %%% hill function stress term
    k1 = 5;
    fsig = (y(6)^k1)/(1 + y(6)^k1);
    
    %%% hill function traction force term
    N0 = 2;
    k2 = 5;
    fn = (y(1)^k2)/(N0^k2 + y(1)^k2);
    
    % y(1) = n, y(2) = n_x, y(3) = rho, y(4) = rho_x, y(5) = u, y(6) =
    % sigma
    f = [y(2)
        y(1)-1-a*fsig
        y(4)
        (1/D2)*(d2*y(3)-y(1))
        (1/E1)*(y(6)-taup2*y(3)*fn)
        0];
end

%% Helper functions

function res = bcfcn(ya,yb) % boundary conditions
res = [yb(1)-1
       ya(2)
       yb(3)-1
       ya(4)
       ya(5)
       yb(5)];
end

function g = lowguess(~)
g = [1
    0
    1
    0
    0
    0];
end

function g = highguess(~)
g = [5
     0
     1
     0
     0
     5];
end