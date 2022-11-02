function [t5pp,t6pp] = SecondOrder_L2(r2,r26,r6,r5,r15,t2,t3,t4,t5,t6,t3p,t4p,t5p,t6p,t3pp,t4pp)
% Second Order Kinematic Coefficient Analysis, Loop 2 - find t5pp and t6pp 
% (theta 5 double prime theta 6 double prime)

%Preallocating Space
t5pp = ones(1,length(t2));
t6pp = ones(1,length(t2));

for i = 1:length(t2)
    
    J = [-r5*sin(t5(i)), -r6*sin(t6(i));
          r5*cos(t5(i)), r6*cos(t6(i))];   
    
    RightSide = [r2*cos(t2(i))+r26*sin(t3(i))*t3pp(i)+r26*cos(t3(i))*t3p(i)^2+r6*cos(t6(i))*t6p(i)^2+...
                 r5*cos(t5(i))*t5p(i)^2+r15*sin(t4(i))*t4pp(i)+r15*cos(t4(i))*t4p(i)^2;
                 r2*sin(t2(i))-r26*cos(t3(i))*t3pp(i)+r26*sin(t3(i))*t3p(i)^2+r6*sin(t6(i))*t6p(i)^2+...
                 r5*sin(t5(i))*t5p(i)^2-r15*cos(t4(i))*t4pp(i)+r15*sin(t4(i))*t4p(i)^2];
    
    Solution = J\RightSide;
    
    t5pp(i) = Solution(1,1);
    t6pp(i) = Solution(2,1);
end

end