function [x,y] = Simulations_1DMix(in_stress, in_traction, in_IC, in_K, in_a, in_taup)
% clc
% close all

% set(0,'DefaultAxesFontName', 'Times New Roman')
% set(0,'DefaultAxesFontSize', 18)
% set(0,'defaultaxeslinewidth',1)
% set(0,'defaultpatchlinewidth',1)
% set(0,'defaultlinelinewidth',2)
% set(0,'defaultTextInterpreter','latex')

%% Numerical set up
par.stress = in_stress;          % 0 for no stress related cell recruitment, 1 for linear, 2 for Hill function
par.traction = in_traction;      % 1 for traction force linear wrt. cells, 2 for Hill function
par.IC = in_IC;                  % 0 for low IC, 1 for high IC
par.K = in_K;                     % Number of spatial grid cells
par.L = 1;                     % Domain length
par.Nright = 1;                  % Right dirichlet BC for cells
par.Pright = 1;                  % Right dirichlet BC for collagen
par.Uleft = 0;                   % Left dirichlet BC for displacement
par.Uright = 0;                  % Right dirichlet BC for displacement
par.x = linspace(0,par.L,par.K); % Discretise spatial domain
par.dx = par.x(2)-par.x(1);      % Cell size
t0 = 0;                          % Initial time
tf = 100;                        % Final time
tspan = linspace(t0,tf,1001);    % Time span

x = linspace(0,par.L,par.K);     % Discretise spatial domain, used for video generation
dt = tspan(2)-tspan(1);          % Time step, used for video generation

%%% For generating the parameter plane
par.a = in_a;
par.taup = in_taup;

%% Initial conditions - eq.(28)
% steadystate = [par.Nright*ones(par.K,1); par.Pright*ones(par.K,1); zeros(par.K,1)];
f = 1 + par.IC*exp(-((par.x)./0.4).^2);
% f1 = 1 + 0.1*exp(-((par.x)./0.4).^2);
n0 = par.Nright.*(f.*ones(1,par.K)).';
p0 = par.Pright.*ones(par.K,1);
% p0 = par.Pright.*(f1.*ones(1,par.K)).';
u0 = zeros(par.K,1);
% u0 = 0.2.*(f1.*(ones(1,par.K))).';
y0 = [n0;p0;u0];

%% Solve with ODE15i
%%% Compute consistent yp0 
% res = @(y,yp)(norm(mechanochemical(y,yp,par)));
% disp(['residuum of steady state = ' ...
%     num2str(res(steadystate,0*steadystate), '%15.10e')]);
opt = odeset('RelTol', 10.0^(-7), 'AbsTol', 10.0^(-7));
nfixed = zeros(par.K,1);
pfixed = zeros(par.K,1);
ufixed = zeros(par.K,1);
[y0,yp0,~] = decic(@(t,y,yp)(mechanochemical(y,yp,par)), t0, ...
    y0, [nfixed; pfixed; ufixed], ...
    [zeros(3*par.K,1)], zeros(3*par.K,1), opt);
% disp(['residuum (from decic) of IC = ' num2str(resnorm, '%15.10e')]);
% disp(['residuum (from res()) of IC = ' num2str(res(y0,yp0), '%15.10e')]);

%% Solve
tic
% [t,y] = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
sol = ode15i(@(t,y,yp)(mechanochemical(y,yp,par)),tspan,y0,yp0);
tfinal = sol.stats.tfinal;
t = linspace(t0,tfinal,tfinal/dt+1);
[y,yp] = deval(sol,t);
toc

%%% Save computed solution to file
filename = ['saved_y1D_' num2str(par.K)];
save(filename, 't', 'y', 'par', 'x', 'yp');

%% Plot
% video_on = true; % Record video: YES (true), NO (false)
% video_filename = [filename '.avi'];
% plot_solution(x,y,yp,t,par,video_on,video_filename);
end

%% Main function implementing the model
function f = mechanochemical(y,yp,par)
    %%% Parameter values - shared with all versions
    D2 = 0.1;
    d2 = 1;

    %%% Parameter values - depends on in_stress
    if par.stress == 0
        E1 = 1;
        eta1 = 1e-8;
    elseif par.stress == 1 || par.stress == 2
        E1 = 100;
        eta1 = 1e-6;
    end
    
    %%% Reshape input vectors
    [n,p,u] = deal(y(1:par.K),y(par.K+1:2*par.K),y(2*par.K+1:3*par.K));
    ntilde = [n(2); n; par.Nright]; 
    ptilde = [p(2); p; par.Pright]; 
    utilde = [0; u; 0];
    [np,pp,up] = deal(yp(1:par.K),yp(par.K+1:2*par.K),...
         yp(2*par.K+1:3*par.K));
    uptilde = [0; up; 0];

    %%% Traction force term
    if par.traction == 1        % Linear traction force
        taup1 = par.taup;
        Tr = taup1.*p.*n;
        Trtilde = [Tr(2); Tr; taup1*ptilde(end)*ntilde(end)];
    elseif par.traction == 2    % Hill function traction force
        taup2 = par.taup;
        N0 = 2;
        k2 = 5;
        hn = (n.^k2)./(N0^k2*ones(size(n)) + n.^k2);
        Tr = taup2.*p.*hn;
        n2 = (ntilde(end)^k2)/(N0^k2 + ntilde(end)^k2);
        Trtilde = [Tr(2); Tr; taup2*ptilde(end)*n2];
    end

    %%% Stress related cell recruitment
    sig = eta1*Mx(uptilde, par) + E1*Mx(utilde,par) + Tr;
    % disp(sig(1));
    if par.stress == 0
        fsig = 0;
    elseif par.stress == 1
        fsig = sig;
    elseif par.stress == 2
        a = par.a;
        k1 = 5;
        fsig = a*(sig.^k1)./(ones(size(sig)) + sig.^k1);
    end
    
    %%% Equation for n
    fn = np - Mxx(ntilde, par) + MA(n,up,par) - fsig - ones(size(n)) + n;

    %%% Equation for p
    fp = pp - D2*Mxx(ptilde,par) + MA(p,up,par) - n + d2*p;

    %%% Equation for u 
    fu = eta1*Mxx(uptilde, par) + E1*Mxx(utilde,par) + Mx(Trtilde, par);

    %%% Full system - eq.(S.1)
    f = [fn; fp; fu];
end

%% Annexed functions

%%% Compute first order derivative on grid cell interfaces 
%%% second order approximation, central
function dx1 = Mx(y,par)
    persistent Mdx;     % Mdx is a K x K+2 matrix
    c1 = [-1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdx = (0.5/par.dx)*[c1, diag(ones(par.K-1,1),1)-diag(ones(par.K-1,1),-1) ,cn];
    dx1 = Mdx*y;
end

%%% Compute second order derivative on grid cell interfaces 
%%% second order approximation, central
function dx2 = Mxx(y,par)
    persistent Mdxx;     % Mdx is a K x K+2 matrix
    c1 = [1; zeros(par.K-1,1)];
    cn = [zeros(par.K-1,1); 1];
    Mdxx = (1/par.dx/par.dx)*[c1, -2*eye(par.K) + diag(ones(par.K-1,1),1) + diag(ones(par.K-1,1),-1) ,cn];
    dx2 = Mdxx*y;
end

%%% Compute flux gradient using first order upwinding without center
%%% averaging
%%% Take material at K grid cells
function fluxdiffx1 = MA(y, vel, par)
    flux = zeros(size(y));

    for i=2:size(flux)
        if vel(i)>0
            flux(i) = y(i-1)*vel(i);
        else
            flux(i) = y(i)*vel(i);
        end
    end

    flux(1) = vel(1)*y(1);

    % compute flux difference per grid cell - def.(S.7) and (S.4)
    fluxdiffx1 = Mx([0; flux; 0],par);
    fluxdiffx1(1) = (flux(2)-flux(1))/par.dx;
    fluxdiffx1(end) = (flux(end)-flux(end-1))/par.dx;

    % vx = Mx([0;vel;0],par);
    % fluxdiffx1(end) = vx(end)*1;

    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is 0
    % fluxdiffx1 = [(1/par.dx)*(flux(2:end)-flux(1:(end-1))); 0];

    % compute flux difference per grid cell via first order FD, assuming
    % flux gradient at the right boundary is the same as to its left 
    % fluxdiffx1 = (1/par.dx)*(flux(2:end)-flux(1:(end-1)));
    % fluxdiffx1 = [fluxdiffx1; fluxdiffx1(end)];
end

%%% Plot solution 
function plot_solution(x,y,yp,t,par,video_on,video_filename)
    if video_on % Initialise video
        vid = VideoWriter(video_filename);
        open(vid);
        figure('Units','normalized','Position',[0 0 0.5 0.45])
    end

    for i=1:length(t)
        clf
        n = [y(1:par.K,i)];
        p = [y(par.K+1:2*par.K,i)];
        u = [y(2*par.K+1:3*par.K,i)];
        v = [yp(2*par.K+1:3*par.K,i)];

        subplot(1,4,1)
        plot(x,n)
        title('$n(t,x)$')
        ylim([0.5,max(1.5,max(n))])
        axis square
        subplot(1,4,2)
        plot(x,p)
        title('$\rho(t,x)$')
        % ylim([10,max(30,max(p))])
        ylim([1,2])
        axis square
        subplot(1,4,3)
        plot(x,u)
        title('$u(t,x)$')
        ylim([-0.1,max(0.1,max(u))])
        axis square
        subplot(1,4,4)
        plot(x,v)
        title('$v(t,x)$')
        ylim([-0.1,max(0.1,max(v))])
        axis square
        a = axes;
        t1 = title([' (t=',num2str(t(i)),')']);
        a.Visible = 'off'; 
        t1.Visible = 'on'; 
        drawnow
        if video_on % Record video
            frame = getframe(gcf);
            size(frame.cdata);
            writeVideo(vid,frame);
            pause(0.1)
        end
    end
    if video_on % Close video
        close(vid)
    end
end