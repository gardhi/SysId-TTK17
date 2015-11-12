
% Simulation_Script
% do not run before estimates from Estimation Script has been calculated

tic

% Choose LOG input
LOG = LOG3;
t1 = 1; % LOG3: 500 ;
t2 = length(LOG.t)-1000; % LOG3: 5000 ;

% Assigning variables
t = LOG.t(t1:t2);
time = length(t);
p_c_ds = LOG.p_c_ds(t1:t2);
q_p = LOG.q_p(t1:t2);
q_bpp = LOG.q_bpp(t1:t2);
z_c = LOG.z_c(t1:t2);

% Parameters for estimation
M = 1e+9;

disp(theta)
disp(C1)
disp(C2)

% Simulation
sim SystemSim

% Verification
toc

%plotting
if 0
    
    figure(1); clf(1)
    plot(q_out.time,q_out.signals.values)
    figure(2); clf(2)
    plot(q_c_out.time,q_c_out.signals.values)
    figure(3); clf(3)
    plot(p_c_out.time,p_c_out.signals.values)
    figure(4); clf(4)
    plot(p_p_out.time,p_p_out.signals.values)
    
    
end