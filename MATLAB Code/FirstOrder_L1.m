function [t3p,t4p] = FirstOrder_L1(r2,r3,r4,t2,t3,t4)
% First Order Kinematic Coefficient Analysis, Loop 1
% Find t3p and t4p (theta 3 prime and theta 4 prime)

%Preallocating Space
t3p = ones(1,length(t2));
t4p = ones(1,length(t2));

for i = 1:length(t2)
    
    J = [-r3*sin(t3(i)), -r4*sin(t4(i));
          r3*cos(t3(i)), r4*cos(t4(i))];
    
    RightSide = [r2*sin(t2(i));-r2*cos(t2(i))];
    
    Solution = J\RightSide;
    
    t3p(i) = Solution(1,1);
    t4p(i) = Solution(2,1);
end

end