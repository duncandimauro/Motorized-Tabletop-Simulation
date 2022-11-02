function [t5,t6] = pos_L2(r2, r26, r6, r5, r15, r1, t2, t3, t4, formation)
%Second Vector Loop Position Analysis - find t5 and t6 for all t2

%Note that "L2" means loop 2

%Preallocating Space
t5 = ones(1, length(t2));
t6 = ones(1, length(t2));

if formation == 1 %without the motor
    t5_0 = 1.21*pi;
    t6_0 = 0;
end

if formation == 2 %with the motor
    t5_0 = 4.30817;
    t6_0 = 0;
end

for i = 1:length(t2)
    %set dt5 and dt6 to randomly large number
    dt5 = 100;
    dt6 = 100;
    D = [dt5;dt6];
    
    %Set up variables for iterations
    n = 0;
    t5_1 = t5_0;
    t6_1 = t6_0;
    
    while norm(D,1)>10^(-6)
        
        n = n + 1;
        %Evaluate f1, f2
        f1 = r2*cos(t2(i)) + r26*cos(t3(i)) + r6*cos(t6_1) + r5*cos(t5_1) + r15*cos(t4(i)) - r1;
        f2 = r2*sin(t2(i)) + r26*sin(t3(i)) + r6*sin(t6_1) + r5*sin(t5_1) + r15*sin(t4(i));
        F = [f1;
             f2];

        %Evaluate Jacobian
        
        J = [-r5*sin(t5_1), -r6*sin(t6_1);
              r5*cos(t5_1), r6*cos(t6_1)];
        
        %return dt5, dt6
        D = -J\F;
        dt5 = D(1);
        dt6 = D(2);
        
        %update t5_1, t6_1
        t5_1 = t5_1 + dt5;
        t6_1 = t6_1 + dt6;
        
    end
    
    %Extract t5_1 and t6_1 for solution
    t5(i) = t5_1;
    t6(i) = t6_1;
    
    %Update guess for next iteration
    t5_0 = t5_1;
    t6_0 = t6_1;
end

end