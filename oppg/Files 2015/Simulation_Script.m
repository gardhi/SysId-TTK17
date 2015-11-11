
% Simulation_Script
% do not run before estimates from Estimation Script has been calculated



% Choose LOG input
LOG = LOG1;

% Assigning variables
time = length(LOG.t);
t = LOG.t;
p_c_ds = LOG.p_c_ds;
q_p = LOG.q_p;
q_bpp = LOG.q_bpp;
p_p = LOG.p_p;
z_c = LOG.z_c;

% Simulation
sim SystemSim

% Verification



