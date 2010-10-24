function result = jumpSolver( theta, sonicMach, alfvenMach, gamma)
% Solves the 3D MHD Jump Conditions for an equilibrium shock wave assuming that the pre-shock region
% conforms to zero z components for both velocity and magnetic field as well as unity values for 
% mass and pressure in the pre-shock region. Gamma has also been set to 5/3 in the equations to 
% simplify the numeric complexity, but that could be changed if necessary. Once the solution of the
% system has been found, the solver looks through the numerical results structure, determines which
% solution is the physical one and returns the data structure for that solution. It will also save
% a mat file with the solution for subsequent loading in the Corrugation shock initializer.
%
%>> theta       angle between the shock normal and v/B in the pre-shock region (degrees)   double
%>> sonicMach   value for the sound mach speed for the pre-shock flow.                     double
%>> alfvenMach  value fot the magnetic (alfven) speed in th pre-shock flow.                double
%<< result      solved physical solution.                                                  struct

    %--- Initialization ---%
    if (nargin < 4 || isempty(gamma));         gamma      = 5/3; end
    if (nargin < 3 || isempty(alfvenMach));    alfvenMach = 1/8; end
    if (nargin < 2 || isempty(sonicMach));     sonicMach = 10;   end
    if (nargin < 1 || isempty(theta));         theta = 10;       end

    %--- Define jump conditions ---%
    %   Knowns:
    %            pre: pr, rho, ms, ma, vz=0, bz=0, tang=tan(theta) gamma=gmult*5/3
    %           post: BX = bx (divergance eqn), VZ=0, BZ=0
    %   Uknowns:
    %            pre: vx, vj, bx, bj
    %           post: RHO, VX, VJ, PR, BJ 

    % Convert number input into strings
    gmult    = 3/5*gamma;
    tanTheta = num2str(tan(deg2rad(theta)),'%0.16f');
    ms       = num2str(sonicMach,'%0.16f');
    ma       = num2str(alfvenMach,'%0.16f');

    syms vx vj bx bj RHO VX VJ PR BJ
    
    massEqn = 'vx = RHO*VX';
    % Full version of the momentum equation. Because bx = BX it can be simplified.
    % momXEqn = 'rho*vx^2 - bx^2 + p + (bx^2+bj^2+bz^2)/2 = RHO*VX^2 - BX^2 + P + (BX^2+BJ^2+BZ^2)/2';
    momXEqn = 'vx^2 + 1 + (bj^2)/2 = RHO*VX^2 + PR + (BJ^2)/2';
    momYEqn = 'vx*vj - bx*bj = RHO*VX*VJ - bx*BJ';
    enerEqn = sprintf(['(%s + (vx^2+vj^2)/2 + bx^2+bj^2)*vx - (bx*vx + bj*vj)*bx = ' ...
                       '(%s*PR + RHO*(VX^2+VJ^2)/2 + bx^2+BJ^2)*VX - (bx*VX + BJ*VJ)*bx'], ...
                        num2str(gmult*5/2, '%0.9g'), num2str(gmult*5/2, '%0.9g'));
    magYEqn = 'vx*bj - bx*vj = VX*BJ - bx*VJ';

    %--- Equations that include numeric input ---%
    sndMachEqn = sprintf('%s = vx*(%s)^(1/2)',ms, num2str(3/(gmult*5), '%0.10g'));
    alfMachEqn = sprintf('%s = vx/bx',ma);
    magEqn     = sprintf('bj/bx = %s',tanTheta);
    velEqn     = sprintf('vj/vx = %s',tanTheta);
    
    %--- Solve and simplify system ---%
    fprintf('Solving symbolic jump system...');
    S = solve(massEqn, momXEqn, momYEqn, enerEqn, magYEqn, magEqn, sndMachEqn, alfMachEqn, velEqn);
    fprintf('complete.\n');

    %--- Print all of the solutions for viewing ---%
    for i=1:length(S.RHO)
        fprintf('\n--- Solution: %d ----------------------------\n',i);
        fprintf('Mass:        %s | %s\n', num2str(1,'%0.3g'),                 num2str(double(S.RHO(i)),'%0.3g') );
        fprintf('Pressure:    %s | %s\n', num2str(1,'%0.3g'),                 num2str(double(S.PR(i)),'%0.3g') );
        fprintf('X Velocity:  %s | %s\n', num2str(double(S.vx(i)),'%0.3g'),   num2str(double(S.VX(i)),'%0.3g') );
        fprintf('Y Velocity:  %s | %s\n', num2str(double(S.vj(i)),'%0.3g'),   num2str(double(S.VJ(i)),'%0.3g') );
        fprintf('X Magnet:    %s | %s\n', num2str(double(S.bx(i)),'%0.3g'),   num2str(double(S.bx(i)),'%0.3g') );
        fprintf('Y Magnet:    %s | %s\n\n', num2str(double(S.bj(i)),'%0.3g'), num2str(double(S.BJ(i)),'%0.3g') );
    end

    %--- Find physical solution ---%
    %        Use the post-shock mass density values to determine which of the solutions is physical. For
    %        a physical solution RHO must be both real and greater than the pre-shock mass density, 1.
    index = 0;
    for i=1:length(S.RHO)
        if (isreal(S.RHO(i)) && double(S.RHO(i)) > 1); index = i; break; end
    end
    
    if index == 0
        fprintf('Unable to find a valid solution for this configuration. Operation aborted\n.');
        result = [];
        return
    end

    fprintf('Found %d to be the physical solution.\n',index);

    result.mass       = [1, double(S.RHO(index))];
    result.pressure   = [1, double(S.PR(index))];
    result.velocity   = [double(S.vx(index)), double(S.VX(index)); double(S.vj(index)), double(S.VJ(index)); 0, 0];
    result.magnet     = [double(S.bx(index)), double(S.bx(index)); double(S.bj(index)), double(S.BJ(index)); 0, 0];
    result.theta      = theta;
    result.sonicMach  = sonicMach;
    result.alfvenMach = alfvenMach;

    gammaParameter    = '';
    if gamma ~= 5/3
        gammaParameter = ['-g', strrep(num2str(gamma, '%0.3g'), '.', '_')];
    end
    
    fileName          = sprintf('ics-t%s-ms%s-ma%s%s.mat', ...
                                strrep( num2str(theta,'%0.3g'), '.', '_'), ...
                                strrep( num2str(sonicMach,'%0.3g'), '.', '_'), ...
                                strrep( num2str(alfvenMach,'%g'), '.', '_'), ...
                                gammaParameter);
           
    save(fileName, 'result');
    fprintf('Solution saved to: %s\n',fileName);
    fprintf('Solver operations complete.\n\n');
end
