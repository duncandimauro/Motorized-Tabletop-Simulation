function [t5p,t6p] = FirstOrder_L2(r2,r26,r6,r5,r15,t2,t3,t4,t5,t6,t3p,t4p)
% First Order Kinematic Coefficient Analysis, Loop 2
% Find t5p and t6p (theta 5 prime and theta 6 prime)

%Preallocating Space
t5p = ones(1,length(t2));
t6p = ones(1,length(t2));

for i = 1:length(t2)
    
    J = [-r5*sin(t5(i)), -r6*sin(t6(i));
          r5*cos(t5(i)), r6*cos(t6(i))];
    
    RightSide = [r2*sin(t2(i))+r26*sin(t3(i))*t3p(i)+r15*sin(t4(i))*t4p(i);
                -r2*cos(t2(i))-r26*cos(t3(i))*t3p(i)-r15*cos(t4(i))*t4p(i)];
    
    Solution = J\RightSide;
    
    t5p(i) = Solution(1,1);
    t6p(i) = Solution(2,1);
end

end