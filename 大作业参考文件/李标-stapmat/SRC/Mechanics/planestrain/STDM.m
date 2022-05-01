function [B,DET]=STDM (R,S,XX,B)
XJ = zeros(2, 2, 'double'); 
RP = 1.0 + R;
SP = 1.0 + S;
RM = 1.0 - R;
SM = 1.0 - S;

P(1,1) = 0.25* SP;
P(1,2) = - P(1,1);
P(1,3) = - 0.25* SM;
P(1,4) = - P(1,3);
P(2,1) = 0.25* RP;
P(2,2) = 0.25* RM;
P(2,3) = - P(2,2);
P(2,4) = - P(2,1);

for I=1:2
    for J=1:2
        DUM = 0.0;
        for K=1:4
            DUM = DUM + P(I,K)*XX(J,K);
        end
        XJ(I,J )= DUM; 
    end      
end
DET = XJ(1,1)* XJ(2,2) - XJ(2,1)* XJ(1,2);
DUM=1.0/DET;
XJI(1,1) =  XJ(2,2)* DUM;
XJI(1,2) = -XJ(1,2)* DUM;
XJI(2,1) = -XJ(2,1)* DUM;
XJI(2,2) =  XJ(1,1)* DUM;
for K=1:4
    K2=2*K;
    B(1,K2-1) = 0;
    B(2,K2  ) = 0;
    for I=1:2
        B(1,K2-1) = B(1,K2-1) + XJI(1,I) * P(I,K);
        B(2,K2  ) = B(2,K2) + XJI(2,I) * P(I,K);
    end
    B(3,K2  ) = B(1,K2-1); 
    B(3,K2-1) = B(2,K2  );
end
end