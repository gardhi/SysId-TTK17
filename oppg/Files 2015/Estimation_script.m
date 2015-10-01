% SysID Computer Exercise

h = 1730; %m; h = h_d = h_a
M = 10e9; %kg/m^4
g = 9.8;

%% finding f_a and rho_a

unNaNed_p_dh3 = LOG3.p_dh(~isnan(LOG3.p_dh));
unNaNed_p_c3 = LOG3.p_c(~isnan(LOG3.p_dh));
unNaNed_t3 = LOG3.t(~isnan(LOG3.p_dh));
p_hdiff = unNaNed_p_dh3-unNaNed_p_c3;
avg_p_hdiff = sum(p_hdiff)/length(p_hdiff);

figure(6); clf(6)
plot(unNaNed_t3,unNaNed_p_dh3); hold on
plot(LOG3.t,LOG3.p_c)
plot(unNaNed_t3,p_hdiff)
line([0 length(LOG3.t)],[avg_p_hdiff avg_p_hdiff])
legend('p_{dh}','p_c','p_{dh}-p_c')

% parameter attained
rho_a = (1/g*h)*avg_p_hdiff;


%% Steady state analysis LOG1

nomZ_c = LOG1.z_c ./ max(LOG1.z_c);
nomQ_c = LOG1.q_c ./ max(LOG1.q_c);
nomQ_p = LOG1.q_p ./ max(LOG1.q_p);
nomP_p = LOG1.p_p ./ max(LOG1.p_p);
nomP_c = LOG1.p_c ./ max(LOG1.p_c);
nomP_c_ds = LOG1.p_c_ds ./ max(LOG1.p_c_ds);
nomP_dh = LOG1.p_dh(~isnan(LOG1.p_dh)) ./ max(LOG1.p_dh);
unNaNed_t1 = LOG1.t(~isnan(LOG1.p_dh));

figure(1); clf(1)
plot(LOG1.t,nomZ_c, '--'); hold on;
plot(LOG1.t,nomQ_c);
plot(LOG1.t,nomQ_p);
plot(LOG1.t,LOG1.q_bpp);
plot(LOG1.t,nomP_p);
plot(LOG1.t,nomP_c);
plot(LOG1.t,nomP_c_ds);
plot(unNaNed_t1, nomP_dh);
title('LOG1 nominal values of system (overview)')
legend('z_c','q_c','q_p','q_{bpp}','p_p','p_c','p_c_{ds}','p_{dh}')


%% Steady state analysis LOG2

nomZ_c = LOG2.z_c / max(LOG2.z_c);
nomQ_c = LOG2.q_c / max(LOG2.q_c);
nomQ_p = LOG2.q_p / max(LOG2.q_p);
nomQ_bpp = LOG2.q_bpp / max(LOG2.q_bpp);
nomP_p = LOG2.p_p / max(LOG2.p_p);
nomP_c = LOG2.p_c ./ max(LOG2.p_c);

figure(1); clf(1)
plot(LOG2.t,nomZ_c); hold on;
plot(LOG2.t,nomQ_c);
plot(LOG2.t,nomQ_p);
plot(LOG2.t,nomQ_bpp);
plot(LOG2.t,nomP_p);
plot(LOG2.t,nomP_c);
title('LOG2 nominal values of flows and valve opening')
legend('z_c','q_c','q_p','q_{bpp}','p_p','p_c')

figure(2); clf(2);
plot(LOG2.t, abs(LOG2.p_c-LOG2.p_p) ./ max(abs(LOG2.p_c-LOG2.p_p))); hold on;
plot(LOG2.t,nomQ_c);
grid on;

%% Steady state analysis LOG3

nomZ_c = LOG3.z_c(1:5000)./max(LOG3.z_c(1:5000));
nomQ_c = LOG3.q_c(1:5000) ./ max(LOG3.q_bpp(1:5000));
nomQ_bpp = LOG3.q_bpp(1:5000) ./ max(LOG3.q_bpp(1:5000));
nomQ_p = LOG3.q_p(1:5000) / max(LOG3.q_p(1:5000));
nomP_p = LOG3.p_p(1:5000) / max(LOG3.p_p(1:5000));
nomP_c = LOG3.p_c(1:5000) ./ max(LOG3.p_c(1:5000));

figure(1); clf(1)
plot(LOG3.t(1:5000),nomZ_c(1:5000)); hold on;
plot(LOG3.t(1:5000),nomQ_c(1:5000));
plot(LOG3.t(1:5000),nomQ_p(1:5000));
plot(LOG3.t(1:5000),nomQ_bpp(1:5000));
plot(LOG3.t(1:5000),nomP_p(1:5000));
plot(LOG3.t(1:5000),nomP_c(1:5000));
title('LOG3 nominal values of flows and valve opening')
legend('z_c','q_c','q_p','q_{bpp}','p_p','p_c')


%% chekcing z and u correlation

% The actuator position zc is controlled by a low-level
% control loop which attempts to track the command input u. In most cases, this actuator
% dynamics is negligible fast, so that we can take zc = u . However, in the data set provided,
% you may notice that this is not always the case.

figure(3); clf(3);
subplot(3,1,1)
plot(LOG1.t, LOG1.z_c); hold on;
plot(LOG1.t, LOG1.u_c)
subplot(3,1,2)
plot(LOG2.t, LOG2.z_c); hold on;
plot(LOG2.t, LOG2.u_c)
subplot(3,1,3)
plot(LOG3.t, LOG3.z_c); hold on;
plot(LOG3.t, LOG3.u_c)

% conclusion LOG3 z_c =~ u_c

%% Steady states in LOG1 :: q_p = q_c = q

ss_interval_1 = 760:930;
ss_interval_2 = 1350:1500;
% 1st interval data point matching
p_dh1_int1 = LOG1.p_dh(ss_interval_1);
unNaNed_p_dh1_int1 = p_dh1_int1(~isnan(p_dh1_int1));
unNaNed_p_c1_int1 = LOG1.p_c(ss_interval_1);
unNaNed_p_c1_int1 = unNaNed_p_c1_int1(~isnan(p_dh1_int1));
% 2nd interval data point matching
p_dh1_int2 = LOG1.p_dh(ss_interval_2);
unNaNed_p_dh1_int2 = p_dh1_int2(~isnan(p_dh1_int2));
unNaNed_p_c1_int2 = LOG1.p_c(ss_interval_2);
unNaNed_p_c1_int2 = unNaNed_p_c1_int2(~isnan(p_dh1_int2));
% time points
unNaNed_t1_int1 = ss_interval_1(~isnan(p_dh1_int1));
unNaNed_t1_int2 = ss_interval_2(~isnan(p_dh1_int2));

avg_q_c_int1 = sum(LOG1.q_c(ss_interval_1)) / length(LOG1.q_c(ss_interval_1));
avg_q_c_int2 = sum(LOG1.q_c(ss_interval_2)) / length(LOG1.q_c(ss_interval_2));


C_a1 = (1/avg_q_c_int1)*(unNaNed_p_c1_int1 - unNaNed_p_dh1_int1 + rho_a*g*h);
C_a2 = (1/avg_q_c_int2)*(unNaNed_p_c1_int2 - unNaNed_p_dh1_int2 + rho_a*g*h);

C_a1 = sum(C_a1)/length(C_a1)
C_a2 = sum(C_a2)/length(C_a2)

C_a = (C_a1 + C_a2) / 2
