function p = kCm(x)
    m=0.00003;
    g=9.81;
    p = [0, -m*sin(x(2))*sin(x(3)), m*cos(x(2))*cos(x(3)), 0, 0, 0;
         m*(cos(x(1))*sin(x(2))*sin(x(3))+sin(x(1))*cos(x(3))),  m*sin(x(1))*cos(x(2))*sin(x(3)),  m*(sin(x(1))*sin(x(2))*cos(x(3))+cos(x(1))*sin(x(3))), 0, 0, 0;
         m*(-sin(x(1))*sin(x(2))*sin(x(3))+cos(x(1))*cos(x(3))),   m*cos(x(1))*cos(x(2))*sin(x(3)),    m*(cos(x(1))*sin(x(2))*cos(x(3))-sin(x(1))*sin(x(3))), 0, 0, 0;
        0, -g*(cos(x(2))), 0, 0, 0, 0;
        g*(cos(x(1))*cos(x(2))),  g*(-sin(x(1))*sin(x(2))), 0, 0 ,0, 0;
        g*(-sin(x(1))*cos(x(2))), g*(-cos(x(1))*sin(x(2))), 0, 0 0, 0;
        0, 0 ,0 ,1, 0, 0;
        0, 0 ,0 ,0, 1, 0;
        0 ,0 ,0 ,0, 0, 1];
end