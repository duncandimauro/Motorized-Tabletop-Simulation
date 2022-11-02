function [t3,t4] = pos_L1(r1, r2, r3, r4, t2, formation)
%First Vector Loop Position Analysis - find t3 and t4 for all t2

%Note that "L1" means loop 1

%Preallocating Space
t3 = ones(1, length(t2));
t4 = ones(1, length(t2));

if formation == 1 %without the motor
    t3_0 = 0.5*pi;
    t4_0 = 1.75*pi;
end

if formation == 2 %with the motor
    t3_0 = 3.7526;
    t4_0 = 5.5278;
end
    
for i = 1:length(t2)
    %set dt3 and dt4 to randomly large number
    dt3 = 100;
    dt4 = 100;
    D = [dt3;dt4];
    
    %Set up variables for iterations
    n = 0;
    t3_1 = t3_0;
    t4_1 = t4_0;
    
    while norm(D,1)>10^(-6)
        
        n = n + 1;
        %Evaluate f1, f2
        f1 = r2*cos(t2(i)) + r3*cos(t3_1) + r4*cos(t4_1) - r1;
        f2 = r2*sin(t2(i)) + r3*sin(t3_1) + r4*sin(t4_1);
        F = [f1;
             f2];

        %Evaluate Jacobian
        
        J = [-r3*sin(t3_1), -r4*sin(t4_1);
              r3*cos(t3_1), r4*cos(t4_1)];
            
        %return dt3, dt3
        D = -J\F;
        dt3 = D(1);
        dt4 = D(2);
        
        %update t3_1, t4_1
        t3_1 = t3_1 + dt3;
        t4_1 = t4_1 + dt4;
        
    end
    
    %Extract t3_1 and t4_1 for solution
    t3(i) = t3_1;
    t4(i) = t4_1;
    
    %Update guess for next iteration
    t3_0 = t3_1;
    t4_0 = t4_1;
end

end