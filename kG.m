function p = kG(x)
    m=0.00003;
    g=9.81;
    p = [m*cos(x(2))*sin(x(3));
         m*(sin(x(1))*sin(x(2))*sin(x(3))-cos(x(1))*cos(x(3)));
         m*(cos(x(1))*sin(x(2))*sin(x(3))-sin(x(1))*cos(x(3)));
        -g*(sin(x(2)));
        g*(sin(x(1))*cos(x(2)));
        g*(cos(x(1))*cos(x(2)));
        x(4);
        x(5);
        x(6)];
end