function [yid, yhat, idMSE,valMSE] = polyARX(m, na, nb, nk, id_u, id_y, val_u, val_y)
%generate a matrix containing the powers for the regressors in each product
%total of numbers which can be represented in base m with na+nb digits
if m > 1
    numb = 0 : m^(na+nb)-1; %generate the numbers in base 10
    numb = dec2base(numb, m); %convert the numbers to base m

    %split the numbers into vectors of digits
    numb2digitvect = zeros(m^(na+nb), na+nb); 
    for i = 1 : m^(na+nb)
        for j = 1 : na+nb
            numb2digitvect(i, j) = str2double(numb(i, j));
        end
    end

    %select only the combinations which have the sum less than or equal to m
    %such that the power degree of the term will not exceed m
    j = 1;
    for i = 1 : m^(na+nb)
        if sum(numb2digitvect(i, :)) <= m
        exponents(j, :) = numb2digitvect(i, :);
        j = j+1;
        end
    end

    %add the cases in which each term is raised to power m
    exponents = [exponents; m*eye(na+nb)]; %matrix with each row representing a possible combination of exponents
elseif m == 1
    exponents = eye(na+nb);
end

% generate the matrix of the regressors
N = length(id_y);
[columns, ~] = size(exponents);
PHI = ones(N, columns);

for i = 1  : N
    for j = 1 : columns
        for k = 1 : na
            index = i-k; %compute the index of delayed output
            if index > 0 
                PHI(i, j) = PHI(i, j) * id_y(index)^exponents(j, k); 
            else
                %if the index is zero or negative the value is considered to
                %be zero and there is no need to compute the other regressors
                %since the product will remain zero
                PHI(i, j) = 0;  
                break
            end
        end
        
        for k = 1 : nb
            index = i-k-nk+1; %compute the index of delayed input
            if index > 0 && PHI(i, j) ~= 0
                
                PHI(i, j) = PHI(i, j) * id_u(index)^exponents(j, na+k);
            else
                %if the index is 0 or negative or the value of the
                %output regressors(y(k-na)) is zero
                %there is no need to compute the other regressors since the
                %value of the producat will remain zero
                PHI(i, j) = 0;
                break;
            end
        end
    end
end

% apply linear regression to compute the parameters
THETA = PHI\id_y;
yid = PHI*THETA;
% validation
% the value of the approximated output will be computed recursively
N = length(val_y);

%an auxiliary vector will be created with the purpose of computing the
%output on the validation data
v = ones(1, columns);
yhat = zeros(N, 1);

for i = 1 : N
    for j = 1 : columns
        for k = 1 : na
            index = i-k;
            if index > 0
                v(1, j) = v(1, j) * yhat(index)^exponents(j, k);
            else
                v(1,j) = 0;
                break
            end
        end
        
        for k = 1 : nb
            index = i-k -nk+1;
            if index > 0 && v(1, j) ~= 0
                v(1, j) = v(1, j) * val_u(index)^exponents(j, na+k);
            else
                v(1, j) = 0;
                break
            end
        end
    end
    yhat(i) = v*THETA;
    % if the approximated equals zero 
    % the value of the validation output is assigned to it
    % consider it as an initial condition
    if yhat(i) == 0
       yhat(i) = val_y(i);
    end
    v = ones(1, columns);
end
% compute MSE values
idMSE = immse(id_y, yid);
valMSE = immse(val_y, yhat);
end

