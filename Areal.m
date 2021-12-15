function p = Areal(x)
    p = [0, 0, 0, 1, tan(x(2))*sin(x(1)), tan(x(2))*cos(x(1));
         0,  0,  0, 0 ,cos(x(1)),  sin(x(1));
         0,   0,  0, 0, sec(x(2))*sin(x(1)), sec(x(2))*cos(x(1));
        0, 0 ,0 ,0,-x(6)/2, 0;
        0, 0 ,0 ,x(6)/2, 0, 0;
        0 ,0 ,0 ,0, 0, 0];
end