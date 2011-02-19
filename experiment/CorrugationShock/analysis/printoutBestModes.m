function printoutBestModes(grorate, phaserate, groR, phaseR, velX)

grorate   = squeeze(grorate);
phaserate = squeeze(phaserate);
groR      = squeeze(groR);
phaseR    = squeeze(phaseR);

dim0 = size(groR);

grorate = grorate(1:(end/2),:);
phaserate = phaserate(1:(end/2),:);
groR = groR(1:(end/2),:);
phaseR = phaseR(1:(end/2),:);

[Kyprime Kzprime] = ndgrid(1:dim0(1), 1:dim0(2));
Kyprime = Kyprime(1:(end/2),:)-1;
Kzprime = circshift(Kzprime(1:(end/2),:) - floor(dim0(2)/2) - 1,[0 -floor(dim0(2)/2)]);

if dim0(1) <= dim0(2)
    Ky = Kyprime * 200 * pi;
    Kz = Kzprime * 200 * pi * dim0(1) / dim0(2);
else
    Kz = Kzprime * 200 * pi;
    Ky = Kyprime * 200 * pi * dim0(2) / dim0(1);
end

Wre_norm = phaserate ./ (sqrt(Ky.^2 + Kz.^2) * velX);
Wim_norm = grorate   ./ (sqrt(Ky.^2 + Kz.^2) * velX);

growRmin = input('Minimum growth rate correlation to require: ');
phaseRmin = input('Minimum phase rate correlation to require: ');

ind1 = find(groR > growRmin);
ind2 = find((groR > growRmin) & (phaseR > phaseRmin));

OUT = fopen('bestmodes.txt','a');

fprintf(OUT,'Preshock V_x0 = %g\n', velX);
fprintf(OUT,'W normalized by 1/(|K| Vin)\n');
fprintf(OUT,'Identified %i modes with growth correlation > %g:\n', numel(ind1), growRmin);
for x = 1:numel(ind1)
    I = ind1(x);
    fprintf(OUT, '\tMode <%-2i %-2i>: Rgr=%4.3g, Rph=%4.3g, Ky=%-8f, Kz=%-8f, What=%-8.6f+%-8.6fi\n',Kyprime(I), Kzprime(I), groR(I), phaseR(I), Ky(I), Kz(I), Wre_norm(I),Wim_norm(I));
end

fprintf(OUT,'Identified %i modes with growth R > %g and omega_re R > %g:\n', numel(ind2), growRmin, phaseRmin);
for x = 1:numel(ind2)
    I = ind2(x);
    fprintf(OUT, '\tMode <%-2i %-2i>: Rgr=%4.3g, Rph=%4.3g, Ky=%-8f, Kz=%-8f, What=%-8.6f+%-8.6fi\n', Kyprime(I), Kzprime(I), groR(I), phaseR(I), Ky(I), Kz(I), Wre_norm(I),Wim_norm(I));
end

fclose(OUT);

end
