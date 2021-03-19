function plotLinFit(X,Y,yfit,P,Powers,Rsq,xlabl,ylabl)
% Plots data and their polynomial fits. Figure title includes polynomial
% coefficients and r-square
%
%__________________________________________________________________________
% Copyright (C) 2021 Laboratory for Neuroimaging Research
% Written by A. Lutti, 2021.
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

base=10;
figure
plot(X,Y,'.')
hold
plot(sort(X),sort(yfit),'m')
Powers=linspace(0,Powers,Powers+1);

ctr=1;
Str1=['(\alpha_' num2str(Powers(ctr))];
if P(ctr)~=0
    RoundingStr=['1e' num2str(-floor(log(abs(P(ctr)))./log(base))+1)];
    DisplayStr1=['1e' num2str(-floor(log(abs(P(ctr)))./log(base)))];
    DisplayStr2=['e' num2str(floor(log(abs(P(ctr)))./log(base)))];
    Str2=['(' num2str(round(P(ctr)*eval(RoundingStr))/eval(RoundingStr)*eval(DisplayStr1)) DisplayStr2];
else
    Str2=['(' num2str(0)];

end
for ctr=2:size(Powers,2)% (size(P,1)-1):-1:1
    Str1=[Str1,[', \alpha_' num2str(Powers(ctr))]];
    if P(ctr)~=0
        RoundingStr=['1e' num2str(-floor(log(abs(P(ctr)))./log(base))+1)];
        DisplayStr1=['1e' num2str(-floor(log(abs(P(ctr)))./log(base)))];
        DisplayStr2=['e' num2str(floor(log(abs(P(ctr)))./log(base)))];
        Str2=[Str2,[', ' num2str(round(P(ctr)*eval(RoundingStr))/eval(RoundingStr)*eval(DisplayStr1)) DisplayStr2]];
    else
        Str2=[Str2,', ' num2str(0)];
    end
end
Str1=[Str1,') = '];Str2=[Str2,')'];
title([Str1 Str2 '; R^2 = ' num2str(round(Rsq*1e2)/1e2)])
ylabel(ylabl);xlabel(xlabl);

end
