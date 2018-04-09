
tspan=[0 10];
x0=zeros(1,6);
x0(1)=pi/20;
x0(2)=-x0(1);
x0(3)=pi/6;
x0(4)=4;
options=odeset('Events',@switch_events);
t=[];
x=[];
step=100;

%x0=[0.3827 0.0242 3.1944 3.2500 3.4397 10.9786];
for i=1:step
    [t_each_step,x_each_step,te,xe,ie]=ode45(@(t,x) dynamics(t,x),tspan,x0,options);
    if isempty(t)
        t=t_each_step;
        x=x_each_step;
    else
        t_each_step=t_each_step+t(end);
        t=[t;t_each_step];
        x=[x;x_each_step];
    end
    [x_after,z_after,~]=impact(x_each_step(end,:));
    x0=x_after(1:6);

    xixi=1;
end

hold on
for i=1:2
    plot(t,x(:,i))
end

%% compute Jacobian, in this case, it is a scalar.
q1_star_minus=1.418944535884725;
perturbation=0.3;
lambda=-1;
poincare=zeros(2,6);
for i=1:2
    x0(1)=0.382699081698739;
    x0(2)=-0.382649160429339;
    x0(3)=0.392781616189557;
    x0(4)=q1_star_minus+(lambda)^i*perturbation;
    x0(5)=-1.418778733622924;
    x0(6)=-3.093160702753254e-04;

    [poincare(i,:),~,delta_qdot]=impact(x0);
end

jacobian=(poincare(2,4)-poincare(1,4))/2*perturbation;

