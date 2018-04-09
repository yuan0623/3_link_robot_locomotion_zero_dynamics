%{
[r,m,Mh,Mt,L,~]=model_params_three_link_numerical_value;

th1=0.1; th2=0.1; th3=0.1;
De=zeros(5,5);
De(1,1)=(r^2*(4*Mh + 4*Mt + 5*m))/4;
De(1,2)=-(m*r^2*(cos(th1)*cos(th2) + sin(th1)*sin(th2)))/2;
De(1,3)=L*Mt*r*(cos(th1)*cos(th3) + sin(th1)*sin(th3));
De(1,4)=(r*cos(th1)*(2*Mh + 2*Mt + 3*m))/2;
De(1,5)=-(r*sin(th1)*(2*Mh + 2*Mt + 3*m))/2;
De(2,1)=-(m*r^2*(cos(th1)*cos(th2) + sin(th1)*sin(th2)))/2;
De(2,2)=(m*r^2)/4;
De(2,4)=-(m*r*cos(th2))/2;
De(2,5)=(m*r*sin(th2))/2;
De(3,1)=L*Mt*r*(cos(th1)*cos(th3) + sin(th1)*sin(th3));
De(3,3)=L^2*Mt;
De(3,4)=L*Mt*cos(th3);
De(3,5)=-L*Mt*sin(th3);
De(4,1)=(r*cos(th1)*(2*Mh + 2*Mt + 3*m))/2;
De(4,2)=-(m*r*cos(th2))/2;
De(4,3)=L*Mt*cos(th3);
De(4,4)=Mh + Mt + 2*m;
De(5,1)=-(r*sin(th1)*(2*Mh + 2*Mt + 3*m))/2;
De(5,2)=(m*r*sin(th2))/2;
De(5,3)=-L*Mt*sin(th3);
De(5,5)=Mh + Mt + 2*m;


theta1=th1; theta2=th2; theta3=th3;
l=L;
De1=zeros(5,5);
De1(1,1)=(5/4*m+Mh+Mt)*r^2;
De1(1,2)=-1/2*m*r^2*cos(theta1-theta2);
De1(1,3)=Mt*r*l*cos(theta1-theta3);
De1(1,4)=(3/2*m+Mh+Mt)*r*cos(theta1);
De1(1,5)=-(3/2*m+Mh+Mt)*r*sin(theta1);
De1(2,1)=De(1,2);
De1(2,2)=1/4*m*r^2;
De1(2,3)=0;
De1(2,4)=-1/2*m*r*cos(theta2);
De1(2,5)=1/2*m*r*sin(theta2);
De1(3,1)=De(1,3);
De1(3,2)=De(2,3);
De1(3,3)=Mt*l^2;
De1(3,4)=Mt*l*cos(theta3);
De1(3,5)=-Mt*l*sin(theta3);
De1(4,1)=De(1,4);
De1(4,2)=De(2,4);
De1(4,3)=De(3,4);
De1(4,4)=2*m+Mh+Mt;
De1(4,5)=0;
De1(5,1)=De(1,5);
De1(5,2)=De(2,5);
De1(5,3)=De(3,5);
De1(5,4)=De(4,5);
De1(5,5)=2*m+Mh+Mt;
%}
syms alpha11 alpha10 alpha21 alpha20 alpha1M alpha1M_minus_1 alpha2M alpha2M_minus_1

H=[0 1 0;0 0 1;1 0 0];
R=[0 1 0;1 0 0;0 0 1];
q1_dot_plus=0.9648;
q1_dot_minus=1.6368;
left_HS=[3.8348*(alpha11-alpha10);3.8348*(alpha21-alpha20);1]*q1_dot_plus;
right_HS=H*delta_qdot/H*[3.8348*(alpha1M-alpha1M_minus_1);3.8348*(alpha2M-alpha2M_minus_1);1]*q1_dot_minus;