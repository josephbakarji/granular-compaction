function [flux, Sm] = Fhllc(Ul, Ur, varnum)
% Compute flux and HLLC godunov-method
% Assumes 7 variables: varnum takes values 1 to 7.
% Ul is value at left grid point
% Ur is value at right grid point.
% Sm is the mean characteristic velocity.
% flux is the resulting flux


global gam1 gam2 pinf1 pinf2 c1 c2

ul      = u(Ul);
rho1l   = rho1(Ul);
rho2l   = rho2(Ul);
El      = E(Ul);

ur      = u(Ur);
rho1r   = rho1(Ur);
rho2r   = rho2(Ur);
Er      = E(Ur);

cl = c(Ul, c1, c2);
cr = c(Ur, c1, c2);
rhol = rho(Ul);
rhor = rho(Ur);

pr = Ur(8);
pl = Ul(8);


Sr = max(ul + cl, ur + cr);
Sl = min(ul - cl, ur - cr);
Sm = ((rhol*ul^2 + pl) - (rhor*ur^2+pr) - Sl*(rhol*ul) + Sr*(rhor*ur))...
    ./(rhol*ul - rhor*ur - Sl*rhol + Sr*rhor);



%% Solver


if(varnum == 1)
    Fl = ul*Ul(1);
    Fr = ur*Ur(1);
    Uls = arhoS(Ul(1), Sl, Sm, ul);
    Urs = arhoS(Ur(1), Sr, Sm, ur);
    flux = hllc(Ul(1), Ur(1), Uls, Urs, Fl, Fr, Sl, Sr, Sm);


elseif(varnum == 2)
    Fl = ul*Ul(2);
    Fr = ur*Ur(2);
    Uls = arhoS(Ul(2), Sl, Sm, ul);
    Urs = arhoS(Ur(2), Sr, Sm, ur);
    flux = hllc(Ul(2), Ur(2), Uls, Urs, Fl, Fr, Sl, Sr, Sm);    


elseif(varnum == 3)

arhoS1r = arhoS(Ur(1), Sr, Sm, ur);
arhoS1l = arhoS(Ul(1), Sl, Sm, ul);
arhoS2r = arhoS(Ur(2), Sr, Sm, ur);
arhoS2l = arhoS(Ul(2), Sl, Sm, ul);

rhoSr = arhoS1r + arhoS2r;
rhoSl = arhoS1l + arhoS2l;

    Fl = ul^2 * rho(Ul)+pl;
    Fr = ur^2 * rho(Ur)+pr;
    Uls = rhoSl * Sm;
    Urs = rhoSr * Sm;
    flux = hllc(Ul(3), Ur(3), Uls, Urs, Fl, Fr, Sl, Sr, Sm);      


elseif(varnum == 4)

arhoS1r = arhoS(Ur(1), Sr, Sm, ur);
arhoS1l = arhoS(Ul(1), Sl, Sm, ul);
arhoS2r = arhoS(Ur(2), Sr, Sm, ur);
arhoS2l = arhoS(Ul(2), Sl, Sm, ul);

rhoSr = arhoS1r + arhoS2r;
rhoSl = arhoS1l + arhoS2l;

pS = pr + rhor*ur*(ur-Sr) - rhoSr*Sm*(Sm-Sr);

ESr = ES(rhor, Er, ur, Sr, pr, rhoSr, Sm, pS);
ESl = ES(rhol, El, ul, Sl, pl, rhoSl, Sm, pS);

    Fl = (El * rho(Ul) + pl) * ul;
    Fr = (Er * rho(Ur) + pr) * ur;
    Uls = ESl * rhoSl;
    Urs = ESr * rhoSr;
    flux = hllc(Ul(4), Ur(4), Uls, Urs, Fl, Fr, Sl, Sr, Sm);


elseif(varnum == 5)
    if(Sm>=0)
        flux = Sm * Ul(5);
    else
        flux = Sm * Ur(5);
    end
    


elseif(varnum == 6)
rhoS1r = arhoS(rho1r, Sr, Sm, ur);
rhoS1l = arhoS(rho1l, Sl, Sm, ul);

p1l = pk(pinf1, gam1, rho1l, e1(Ul));
p1r = pk(pinf1, gam1, rho1r, e1(Ur));

pS1l = pSk(gam1, p1l, pinf1, rhoS1l, rho1l);
pS1r = pSk(gam1, p1r, pinf1, rhoS1r, rho1r);

eS1r = eSS(gam1, pinf1, pS1r, rhoS1r); 
eS1l = eSS(gam1, pinf1, pS1l, rhoS1l); 

    Fl = ul*Ul(6);
    Fr = ur*Ur(6);
    Uls = arhoS(Ul(1), Sl, Sm, ul) * eS1l;
    Urs = arhoS(Ur(1), Sr, Sm, ur) * eS1r;
    flux = hllc(Ul(6), Ur(6), Uls, Urs, Fl, Fr, Sl, Sr, Sm);
    


elseif(varnum == 7)
rhoS2r = arhoS(rho2r, Sr, Sm, ur);
rhoS2l = arhoS(rho2l, Sl, Sm, ul);

p2l = pk(pinf2, gam2, rho2l, e2(Ul));
p2r = pk(pinf2, gam2, rho2r, e2(Ur));

pS2l = pSk(gam2, p2l, pinf2, rhoS2l, rho2l);
pS2r = pSk(gam2, p2r, pinf2, rhoS2r, rho2r); % BUUUUUG!

eS2r = eSS(gam2, pinf2, pS2r, rhoS2r); 
eS2l = eSS(gam2, pinf2, pS2l, rhoS2l); 

    Fl = ul*Ul(7);
    Fr = ur*Ur(7);
    Uls = arhoS(Ul(2), Sl, Sm, ul) * eS2l;
    Urs = arhoS(Ur(2), Sr, Sm, ur) * eS2r;
    flux = hllc(Ul(7), Ur(7), Uls, Urs, Fl, Fr, Sl, Sr, Sm);
else
    flux = nan;
    disp('incorrect variable input')
    
end

