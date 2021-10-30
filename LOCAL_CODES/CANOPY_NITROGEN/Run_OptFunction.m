%% Canopy Nitrogen distribution

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%
%-------------------------------------------------------------------------%
%   Created by   : Rohit Nandan                                           %
%   Date         : April, 2021                                      %
%-------------------------------------------------------------------------%

% load ('temp_variable_N_temp.mat', 'PARAMS', 'VARIABLES', 'VERTSTRUC');

% K_ec, K_rub_leaf, K_lhc_leaf

lb(1)=MinKec; lb(2)=MinFrub; lb(3)=MinFth;
ub(1)=MaxKec; ub(2)= MaxFrub; ub(3)= MaxFth;

%% Define name of function
fitnessfcn = @CAN_NIT_function;

nvars = 3;

% rng(1, 'twister')
rng default         % For reproducibility

options = optimoptions('gamultiobj','UseParallel', SWITCHES.ParOnV, 'UseVectorized', false, 'PopulationSize', 500);

[x,fval,exitflag,output,population,scores] = gamultiobj(fitnessfcn,nvars, ...
    [],[],[],[],lb,ub,options);

save Optimization_output.mat x fval exitflag output population scores;

save ([workpath, '/LOCAL_CODES/CANOPY_NITROGEN/AuxFiles/NOptS', num2str(si),'Y',...
    num2str(Run_years(yy)), 'D',...
    num2str(fix(tt/24)), '.mat'], 'x', 'fval', 'scores')

% Sum_An = sum(fval,2);
% [Value, OrderNo] = find(Min_An);

Min_An=[];

[c, r]=size(fval);
[cx, rx]=size(x);
% disp(c)

for ii=1:c
    Min_An(ii)=min([fval(ii,1), fval(ii,2)]);
end

[Value, OrderNo] = max(Min_An);

% disp(OrderNo)

if c==cx && c>=OrderNo
    if isempty(OrderNo)
        K_ec=0.8;
        K_rub_leaf=0.13;
        K_lhc_leaf=0.25;
    else
        K_ec=x(OrderNo , 1);
        K_rub_leaf=x(OrderNo , 2);
        K_lhc_leaf=x(OrderNo , 3);
    end
    
    disp('done')
    
else
    K_ec=x(end , 1);
    K_rub_leaf=x(end , 2);
    K_lhc_leaf=x(end , 3);
end

save ('NCAN_opt_variables.mat', 'K_ec', 'K_rub_leaf', 'K_lhc_leaf');


