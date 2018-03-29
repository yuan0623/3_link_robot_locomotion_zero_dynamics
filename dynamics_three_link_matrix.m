function [D,C,G,B]=dynamics_three_link_matrix(x)


[r,m,Mh,Mt,L,g]=model_params_three_link_numerical_value;

th1=x(1); th2=x(2); th3=x(3);
dth1=x(4); dth2=x(5); dth3=x(6);

% D matrix
D=zeros(3);
D(1,1)=Mh*r^2 + Mt*r^2 + (5*m*r^2)/4;
D(1,2)=-(m*r^2*cos(th1 - th2))/2;
D(1,3)=L*Mt*r*cos(th1 - th3);
D(2,1)=D(1,2);
D(2,2)=(m*r^2)/4;
D(3,1)=D(1,3);
D(3,3)=L^2*Mt;

% C matrix
C=zeros(3);
C(1,2)=-(dth2*m*r^2*sin(th1 - th2))/2;
C(1,3)=L*Mt*dth3*r*sin(th1 - th3);
C(2,1)=(dth1*m*r^2*sin(th1 - th2))/2;
C(3,1)=-L*Mt*dth1*r*sin(th1 - th3);

% G matrix
G=zeros(3,1);
G(1)=- Mh*g*r*sin(th1) - Mt*g*r*sin(th1) - (3*g*m*r*sin(th1))/2;
G(2)=(g*m*r*sin(th2))/2;
G(3)=-L*Mt*g*sin(th3);

% B matrix
B=zeros(3,2);
B(1,1)=-1;
B(2,2)=-1;
B(3,1)=1;
B(3,2)=1;

end