function u=dim_sum_feedback_linearization(x,Fx,Gx)
    th1=x(1);
    th2=x(2);
    th3=x(3);
    dth1=x(4);
    dth2=x(5);
    dth3=x(6);
    %v=[0;0];
    [H,LfH,~]=feedback_linearization(x,Fx,Gx);
    % Bernstein-Bhat controller (uses feedback linearization)
    v=Bhat_and_Bernstein(H,LfH);
    v=v';
    %Kp=1;
    %Kd=0.2;
    %v=-Kp*[th2+th1;th3-pi/6]-Kd*[dth2+dth1;dth3];
    %v=[0;0];
    [th3desired,th1desired,alpha,epsilon]=control_params;
    [D,C,G,B]=dynamics_three_link_matrix(x);
    c_vector=C*[dth1;dth2;dth3]+G;
    dhdq=[0 0 1;
          1 1 0];
    bqqdot=[0;0];
    u=(dhdq/D*B)\(v-bqqdot+dhdq/D*c_vector);
    
    
end