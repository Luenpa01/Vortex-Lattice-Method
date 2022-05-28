function C = Ctrap(y,Ala)
% Distribución de cuerdas pAla.Parametro.ARa ala trapecial

Ala.Parametro.b = (Ala.Parametro.AR*Ala.Parametro.Sw)^0.5;
Cr = 2/(Ala.Parametro.TR+1)*Ala.Parametro.Sw/Ala.Parametro.b;
C = 2*Cr/Ala.Parametro.b*(Ala.Parametro.TR-1)*abs(y)+Cr;

end