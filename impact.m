function [x_after,z_after]=impact(x)
    [r,m,Mh,Mt,l,~]=model_params_three_link_numerical_value;
    
    
    De=zeros(5);
    theta1=x(1);
    theta2=x(2);
    theta3=x(3);
    omega1=x(4);
    omega2=x(5);
    omega3=x(6);
    De(1,1)=(5/4*m+Mh+Mt)*r^2;
    De(1,2)=-1/2*m*r^2*cos(theta1-theta2);
    De(1,3)=Mt*r*l*cos(theta1-theta3);
    De(1,4)=(3/2*m+Mh+Mt)*r*cos(theta1);
    De(1,5)=-(3/2*m+Mh+Mt)*r*sin(theta1);
    De(2,1)=De(1,2);
    De(2,2)=1/4*m*r^2;
    De(2,3)=0;
    De(2,4)=-1/2*m*r*cos(theta2);
    De(2,5)=1/2*m*r*sin(theta2);
    De(3,1)=De(1,3);
    De(3,2)=De(2,3);
    De(3,3)=Mt*l^2;
    De(3,4)=Mt*l*cos(theta3);
    De(3,5)=-Mt*l*sin(theta3);
    De(4,1)=De(1,4);
    De(4,2)=De(2,4);
    De(4,3)=De(3,4);
    De(4,4)=2*m+Mh+Mt;
    De(4,5)=0;
    De(5,1)=De(1,5);
    De(5,2)=De(2,5);
    De(5,3)=De(3,5);
    De(5,4)=De(4,5);
    De(5,5)=2*m+Mh+Mt;
    
    E=zeros(2,5);
    E(1,1)=r*cos(theta1);
    E(1,2)=-r*cos(theta2);
    E(1,4)=1;
    E(2,1)=-r*sin(theta1);
    E(2,2)=r*sin(theta2);
    E(2,5)=1;
    dqe_minus_1=[omega1;omega2;omega3;0;0];
    dqe_plus_1=([De -E';E zeros(2)])\[De*dqe_minus_1;0;0];
    F=dqe_plus_1(6:7);
    dqe_plus_1=dqe_plus_1(1:5);
    x_after=zeros(6,1);
    x_after(1)=theta2;
    x_after(2)=theta1;
    x_after(3)=theta3;
    x_after(4)=dqe_plus_1(2);
    x_after(5)=dqe_plus_1(1);
    x_after(6)=dqe_plus_1(3);
    z_after=dqe_plus_1(4:5);
    
    
    
end