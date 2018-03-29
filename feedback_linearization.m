function [H,LfH,dLfH]=feedback_linearization(x,Fx,Gx)
    a = [0.512 0.073 0.035 -0.819 -2.27 3.26 3.11 1.89];
    th1=x(1);
    th2=x(2);
    th3=x(3);
    dth1=x(4);
    dth2=x(5);
    dth3=x(6);
    a01=a(1); a11=a(2); a21=a(3); a31=a(4);
    a02=a(5); a12=a(6); a22=a(7); a32=a(8);
    
    th1d=pi/8;
    H(1,1)=th3 - a01 - a11*th1 - a21*th1^2 - a31*th1^3;
    H(2,1)=th1 + th2 - (th1 + th1d)*(th1 - th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3);
    
    LfH=zeros(2,1);
    LfH(1,1)=dth3 - dth1*(a11 + 2*a21*th1 + 3*a31*th1^2);
    LfH(2,1)=dth2 - dth1*((th1 - th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3) + (th1 + th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3) + (th1 + th1d)*(th1 - th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) - 1);
    
    
    dLfH=zeros(2,6);
    dLfH(1,1)=-dth1*(2*a21 + 6*a31*th1);
    dLfH(1,4)=- a11 - 2*a21*th1 - 3*a31*th1^2;
    dLfH(1,6)=1;
    dLfH(2,1)=-dth1*(2*a02 + 2*(th1 + th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) + 2*a12*th1 + 2*(th1 - th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) + 2*a22*th1^2 + 2*a32*th1^3 + (th1 + th1d)*(2*a22 + 6*a32*th1)*(th1 - th1d));
    dLfH(2,4)=1 - (th1 + th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3) - (th1 + th1d)*(th1 - th1d)*(a12 + 2*a22*th1 + 3*a32*th1^2) - (th1 - th1d)*(a02 + a12*th1 + a22*th1^2 + a32*th1^3);
    dLfH(2,5)=1;

    %theta3_desired=pi/6;
    %H=[theta3-theta3_desired;theta1+theta2];
    %LfH=[0 0 1 0 0 0;1 1 0 0 0 0]*Fx;
    %LfLfh=0;
end