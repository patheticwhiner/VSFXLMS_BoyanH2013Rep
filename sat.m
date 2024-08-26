function x = sat(x,lower,upper)
    x(x<lower) = lower;
    x(x>upper) = upper;
end