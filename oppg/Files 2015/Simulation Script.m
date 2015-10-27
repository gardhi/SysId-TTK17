

% Simulation Script

% Choose LOG input
LOG = LOG1;
time = length(LOG.t);
t = LOG.t;

q_c = LOG.q_c;
q_bpp = LOG.q_bpp;
p_p = LOG.p_p;

% Get Estimates
Estimation_script;
close all

sim SystemSim



