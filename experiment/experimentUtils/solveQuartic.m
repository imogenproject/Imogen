function result = solveQuartic(a, b, c, d, e)

alpha = -3*b^2/(8*a^2) + c/a;
beta  = .125*b^3/a^3 - .5*b*c/a^2 + d/a;
gamma = -3*b^4/(256*a^4) + c*b^2/(16*a^3) - .25*b*d/a^2 + e/a;

if beta == 0
    % lol not likely
    result(1) = -.25*b/a + sqrt(.5*(-alpha + sqrt(alpha^2 - 4*gamma)));
    result(2) = -.25*b/a + sqrt(.5*(-alpha - sqrt(alpha^2 - 4*gamma)));
    result(3) = -.25*b/a - sqrt(.5*(-alpha + sqrt(alpha^2 - 4*gamma)));
    result(4) = -.25*b/a - sqrt(.5*(-alpha - sqrt(alpha^2 - 4*gamma)));
else
    P = -alpha^2/12 - gamma;
    Q = -(alpha^3) / 108 + alpha*gamma/3 - .125*beta^2;
    R = -Q/2 + sqrt(.25*Q^2 + (P^3)/27);
    U = R^(1/3);
%R alpha U P Q]'

%25*Q^2 + (P^3)/27

    if U ~= 0
        y = -5*alpha/6 + U - P/(3*U);
    else
        y = -5*alpha/6 + U - Q^(1/3);
    end

    W = sqrt(alpha + 2*y);

    result(1) = -.25*b/a + (W - sqrt(-(3*alpha+2*y + 2*beta/W))) / 2;
    result(2) = -.25*b/a + (W + sqrt(-(3*alpha+2*y + 2*beta/W))) / 2;
    result(3) = -.25*b/a - (W - sqrt(-(3*alpha+2*y - 2*beta/W))) / 2;
    result(4) = -.25*b/a - (W + sqrt(-(3*alpha+2*y - 2*beta/W))) / 2;

end

end
