function dx = dynamics(t,x)


%size(x)
[D,C,G,B] = dynamics_three_link_matrix(x);
Fx = [x(4:6);(D)\(-C*x(4:6)-G)];
Gx = [zeros(3,2);(D)\B];


%[H,LfH,dLfH]=feedback_linearization(x,Fx,Gx);
% Bernstein-Bhat controller (uses feedback linearization)
%v=Bhat_and_Bernstein(H,LfH);
%v=v';
% Used for controller that use feedback linearization
%u = (dLfH*Gx)\(v-dLfH*Fx);
u=dim_sum_feedback_linearization(x,Fx,Gx);

dx= Fx+Gx*u;



end