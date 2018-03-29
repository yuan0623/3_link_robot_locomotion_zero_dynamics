function v=Bhat_and_Bernstein(x1,x2)
    alpha=0.9;
    epsilon=0.1;
    phi(1)=x1(1)+(1/2-alpha)*sign(x2(1))*abs(x2(1))^(2-alpha);
    phi(2)=x1(2)+(1/2-alpha)*sign(x2(2))*abs(x2(2))^(2-alpha);
    
    psi(1)=-sign(x2(1))*abs(x2(1))^alpha-sign(phi(1))*abs(phi(1))^(alpha/(2-alpha));
    psi(2)=-sign(x2(2))*abs(x2(2))^alpha-sign(phi(2))*abs(phi(2))^(alpha/(2-alpha));
    v=1/epsilon^3*psi;

end