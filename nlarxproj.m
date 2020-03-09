%%
clear 
close all

%% load data as object of type iddata from the System Identification toolbox
load('_____');

%% set the parameters according to the matlab narx model
m = 2;
na = 1;
nb = 4;
nk = 2;

%%
[yid, yhat, idMSE,valMSE] = polyARX(m, na, nb, nk, id.u, id.y, val.u, val.y);

%%
subplot(211)
plot(id.y);
hold on
plot(yid);
hold off
title(['ID data fit; MSE = ',num2str(idMSE)]);
grid
subplot(212)
plot(val.y);
hold on
plot(yhat);
hold off
title(['VAL data fit; MSE = ',num2str(valMSE)]);
grid
