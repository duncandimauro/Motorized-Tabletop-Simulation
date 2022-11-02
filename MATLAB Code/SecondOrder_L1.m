function [t3pp,t4pp] = SecondOrder_L1(r2,r3,r4,t2,t3,t4,t3p,t4p)
% Second Order Kinematic Coefficient Analysis, Loop 1 - find t3pp and t4pp 
% (theta 3 double prime theta 4 double prime)

%Preallocating Space
t3pp = ones(1,length(t2));
t4pp = ones(1,length(t2));

for i = 1:length(t2)
    
    J = [-r3*sin(t3(i)), -r4*sin(t4(i));
          r3*cos(t3(i)), r4*cos(t4(i))];
    
    RightSide = [r2*cos(t2(i))+r3*cos(t3(i))*t3p(i)^2+r4*cos(t4(i))*t4p(i)^2;
                 r2*sin(t2(i))+r3*sin(t3(i))*t3p(i)^2+r4*sin(t4(i))*t4p(i)^2];
    
    Solution = J\RightSide;
    
    t3pp(i) = Solution(1,1);
    t4pp(i) = Solution(2,1);
end

end