function [value,isterminal,direction] = switch_events(t,x)


    [~,th1desired,~,~] = control_params;
    tolerence1=(th1desired-x(1));
    value(1) = double(~(tolerence1<0.01));         % when stance leg reach angle of th1d    
    isterminal=1;
 
    direction=-1;
end