
format
% SysID Computer Exercise

h = 1730; %m; h = h_d = h_a
M = 10e9; %kg/m^4
g = 9.8;

time_span = (length(LOG3.t)-1)/60

%% finding rho_a

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
rho_a = (1/(g*h))*avg_p_hdiff;


%% Steady state analysis LOG1

nomZ_c = LOG1.z_c ./ max(abs(LOG1.z_c));
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
%plot(LOG2.t, abs(LOG2.p_c-LOG2.p_p) ./ max(abs(LOG2.p_c-LOG2.p_p))); hold on;
plot(LOG2.t,LOG2.q_c);
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





%% Steady State Curve fitting version 2

t1 = 1;
t2 = length(LOG2.t);
t = LOG2.t(t1:t2);

% Moving average smooths out graphs
a = 10;
B = 1/a*ones(a,1);
q = filter(B,1,LOG2.q_c(t1:t2));
p_c = filter(B,1,LOG2.p_c(t1:t2));
p_p = filter(B,1,LOG2.p_p(t1:t2));

figure(3); clf(3)
plot(t, q); hold on
plot(t, p_p)
plot(t, p_c)


% finding saddle points
derivatives = zeros(1,t2);
for i = a:t2-1
    derivatives(i) = q(i+1)-q(i);
end
tollerance = 0.1;
saddle_points = find(abs(derivatives) < tollerance);
t_ss = t(saddle_points);
t_ss = t_ss(10:end);
q_ss = q(saddle_points);
q_ss = q_ss(10:end);
p_c_ss = p_c(saddle_points);
p_c_ss = p_c_ss(10:end);
p_p_ss = p_p(saddle_points);
p_p_ss = p_p_ss(10:end);

plot(t_ss, q_ss, 'ored')

figure(4);clf(4)
plot(q_ss, p_p_ss - p_c_ss - (rho_a*g*h))
hold on
poly = polyfit(q_ss, p_p_ss - p_c_ss - (rho_a*g*h),2);
plot(q_ss, poly(1)*q_ss.^2 + poly(2)*q_ss + poly(3),'ored')

rho_d = -poly(3)/(g*h)
C_a = poly(2)
D_d = poly(1)


%% Estimating g_c
LOG = LOG3;

figure(10); clf(10)
plot(LOG.t, LOG.p_c-LOG.p_c_ds)
t1 = 500; t2 = 5000;

g_c = LOG.q_c(t1:t2) ./ sqrt(LOG.p_c(t1:t2) -LOG.p_c_ds(t1:t2));

figure(12);clf(12)
plot(LOG.z_c(t1:t2), g_c)
hold on

theta = polyfit(LOG.z_c(t1:t2), g_c, 2);
plot(LOG.z_c(t1:t2), theta(1)*LOG.z_c(t1:t2).^2 + theta(2)*LOG.z_c(t1:t2) +theta(3), 'o')

%% Estimating Va/Beta_a version 2

t1 = 800; t2 = 1800;

LOG = LOG3;
p_c = LOG.p_c(t1:t2);
q_bpp = LOG.q_bpp(t1:t2);
q_c = LOG.q_c(t1:t2);
t = LOG.t(t1:t2);

dp_c_dt = diff(p_c);
cut = 200;

figure(3); clf(3)
plot(dp_c_dt(1:end-cut+1), q_bpp(1:end-cut) - q_c(1:end-cut),'.'); hold on
C1_poly = polyfit(dp_c_dt(1:end-cut+1), q_bpp(1:end-cut) - q_c(1:end-cut),1);
plot(dp_c_dt(1:end-cut+1), C1_poly(1)*dp_c_dt(1:end-cut+1)+ C1_poly(2), 'ored')

C1 = C1_poly(1);


%% Estimating V_d/Beta_d version 2

LOG = LOG1;
% LOG2 intervals
%t1 = 500; t2 = 3000;
% LOG1 intervals
t1 = 400; t2 = 600;

p_c = LOG.p_c(t1:t2); 
p_p = LOG.p_p(t1:t2);
q_bpp = LOG.q_bpp(t1:t2);
q_c = LOG.q_c(t1:t2);
q_p = LOG.q_p(t1:t2);
t = LOG.t(t1:t2);

dp_p_dt = diff(p_p);
dp_c_dt = diff(p_c);

C2_poly = polyfit(dp_p_dt, q_p(1:end-1)-(C1*dp_c_dt)+q_bpp(1:end-1)...
    -q_c(1:end-1),1);
figure(5); clf(5)
plot(dp_p_dt, q_p(1:end-1)-(C1*dp_c_dt)+q_bpp(1:end-1)- q_c(1:end-1),'.');
hold on
plot(dp_p_dt, C2_poly(1)*dp_p_dt + C2_poly(2), 'ored')

C2 = C2_poly(1);
%%



% 
% %% Attempt at finding C_a, D_d and rho_d :: q_p = q_c = q
% % Wrong attempt
% 
% ss_interval_1 = 760:930;
% ss_interval_2 = 1350:1500;
% % 1st interval data point matching
% p_dh1_int1 = LOG1.p_dh(ss_interval_1);
% unNaNed_p_dh1_int1 = p_dh1_int1(~isnan(p_dh1_int1));
% unNaNed_p_c1_int1 = LOG1.p_c(ss_interval_1);
% unNaNed_p_c1_int1 = unNaNed_p_c1_int1(~isnan(p_dh1_int1));
% % 2nd interval data point matching
% p_dh1_int2 = LOG1.p_dh(ss_interval_2);
% unNaNed_p_dh1_int2 = p_dh1_int2(~isnan(p_dh1_int2));
% unNaNed_p_c1_int2 = LOG1.p_c(ss_interval_2);
% unNaNed_p_c1_int2 = unNaNed_p_c1_int2(~isnan(p_dh1_int2));
% % time points
% unNaNed_t1_int1 = ss_interval_1(~isnan(p_dh1_int1));
% unNaNed_t1_int2 = ss_interval_2(~isnan(p_dh1_int2));
% 
% avg_q_c_int1 = sum(LOG1.q_c(ss_interval_1)) / length(LOG1.q_c(ss_interval_1));
% avg_q_c_int2 = sum(LOG1.q_c(ss_interval_2)) / length(LOG1.q_c(ss_interval_2));
% 
% 
% C_a1 = (1/avg_q_c_int1)*(unNaNed_p_c1_int1 - unNaNed_p_dh1_int1 + rho_a*g*h);
% C_a2 = (1/avg_q_c_int2)*(unNaNed_p_c1_int2 - unNaNed_p_dh1_int2 + rho_a*g*h);
% 
% C_a1 = sum(C_a1)/length(C_a1);
% C_a2 = sum(C_a2)/length(C_a2);
% 
% C_a = (C_a1 + C_a2) / 2;
% 
% %% Steady state curve fitting version 1.
% % For finding C_a, D_d and rho_d
% 
% q_set = zeros(400,200); % wont need more than one q pr every second point
% p_c_set = zeros(400,200);
% p_p_set = zeros(400,200);
% avg_q_set = zeros(1,length(LOG2.t/2));
% m = 1;
% n = 1;
% set_ongoing = false;
% 
% for i = 10:(length(LOG2.t)-10);
%     
%     if(( abs(LOG2.q_c(i-7)-LOG2.q_c(i)) < 10 ) ...         % checking k1 samples to the left
%             && ( abs(LOG2.q_c(i+7)-LOG2.q_c(i)) < 10 ) )   % checking k2 samples to the right
%         q_set(n,m) = LOG2.q_c(i);                          % if not too different assume steady state
%         p_c_set(n,m) = LOG2.p_c(i);
%         p_p_set(n,m) = LOG2.p_p(i);
%         n = n + 1;
%         set_ongoing = true;
%     else
%         if set_ongoing
%             m = m+1;
%             n = 1;
%             set_ongoing = false;
%         end
%         
%     end
% end
% 
% measurments_for_use = [];
% k = 1;
% for i = 1: length(q_set(1,:))
%     temp = length(find(q_set(:,i)));
%     if temp > 5;
%         measurments_for_use(k) = i;
%         k = k+1;
%     end
% end
% 
% avg_q_measurments = zeros(1,length(measurments_for_use));
% avg_p_c_measurments = zeros(1,length(measurments_for_use));
% avg_p_p_measurments = zeros(1,length(measurments_for_use));
% 
% for i = 1:length(measurments_for_use)
%     this_q_set = q_set(find(q_set(:,measurments_for_use(i))), measurments_for_use(i));
%     avg_q_measurments(i) = sum(this_q_set)/length(this_q_set);
%     this_p_c_set = p_c_set(find(p_c_set(:,measurments_for_use(i))), measurments_for_use(i));
%     avg_p_c_measurments(i) = sum(this_p_c_set)/length(this_p_c_set);
%     this_p_p_set = p_p_set(find(p_p_set(:,measurments_for_use(i))), measurments_for_use(i));
%     avg_p_p_measurments(i) = sum(this_p_p_set)/length(this_p_p_set);
% end   
% 
% figure(11); clf(11);
% plot(avg_q_measurments ,  - (rho_a*g*h) - avg_p_c_measurments + avg_p_p_measurments)
% hold on;
% poly = polyfit(avg_q_measurments , - (rho_a*g*h) - avg_p_c_measurments + avg_p_p_measurments,2);
% plot(avg_q_measurments, poly(1)*avg_q_measurments.^2 + poly(2)*avg_q_measurments + poly(3), 'o')
% hold off;
% % 
% % rho_d = poly(3)/(g*h)
% % C_a = poly(2)
% % D_d = poly(1)
% 


%%

% %% Phi_c(u) and theta_c as a part of q_c(p_c,u_c) and g_c(u)
% % making p_c our y, and u_c our u.....
% 
% % dp_c/dt euler
% for i = 1:(length(LOG3.p_c)-1)
%     p_c_dot = (LOG3.p_c(i+1)-LOG3.p_c) ./ (LOG3.t(i+1) -LOG3.t(i));
% end
% V_a_over_Beta_a = (LOG3.q_bpp - LOG3.q_c)./p_c_dot;
% 
% % Moving average
% a = 200;
% B = 1/a*ones(a,1);
% smoothed_q_bqq = filter(B,1,LOG3.q_bpp);
% smoothed_q_c = filter(B,1, LOG3.q_c);
% 
% figure(1); clf(1)
% plot(LOG3.t, LOG3.q_c); hold on;
% plot(LOG3.t, smoothed_q_c)
% 
% V_a_over_Beta_a_smoothed = (smoothed_q_bqq - smoothed_q_c) ./ p_c_dot;
% 
% figure(2); clf(2)
% plot(LOG3.t, V_a_over_Beta_a_smoothed); hold on;
% plot(LOG3.t, V_a_over_Beta_a)

%%

% %% chekcing z and u correlation
% 
% % The actuator position zc is controlled by a low-level
% % control loop which attempts to track the command input u. In most cases, this actuator
% % dynamics is negligible fast, so that we can take zc = u . However, in the data set provided,
% % you may notice that this is not always the case.
% 
% figure(3); clf(3);
% subplot(3,1,1)
% plot(LOG1.t, LOG1.z_c); hold on;
% plot(LOG1.t, LOG1.u_c)
% subplot(3,1,2)
% plot(LOG2.t, LOG2.z_c); hold on;
% plot(LOG2.t, LOG2.u_c)
% subplot(3,1,3)
% plot(LOG3.t, LOG3.z_c); hold on;
% plot(LOG3.t, LOG3.u_c)
% 
% % conclusion LOG3 z_c =~ u_c


%%
% %% Estimating Va/Beta_a version 1
% 
% t1 = 800; t2 = 1800;
% 
% LOG = LOG3;
% p_c = LOG.p_c(t1:t2);
% q_bpp = LOG.q_bpp(t1:t2);
% q_c = LOG.q_c(t1:t2);
% t = LOG.t(t1:t2);
% 
% figure(4);clf(4)
% q_bpp = smooth(q_bpp,40);
% q_c = smooth(q_c,1000);
% p_c = smooth(p_c,500);
% 
% plot(t,q_bpp); hold on
% plot(t,q_c)
% plot(t,p_c)
% 
% dp_c_dt = diff(p_c);
% % smoothing of derivative
% %  figure(5); clf(5)
% %  plot(t(1:end-1),dp_c_dt); hold on
% % dp_c_dt = smooth(dp_c_dt,50);
% % plot(t(1:end-1),dp_c_dt)
% 
% cut = 200;
% 
% C1_set = (1./dp_c_dt(1:end-cut+1)) .* ( q_bpp(1:end-cut) - q_c(1:end-cut) );
% 
% t = 1:801;
% % % Plotting
% figure(2);clf(2)
% plot(t, C1_set)
% 
% % Avereging out noise
% C1 = sum(C1_set(500:800))/(800-500);
%  



