%* *****************************************************************
%* - Function of STAPMAT in Solver/Dynamics phase                           *
%*                                                                 *
%* - Purpose:                                                      *
%*     求解时间积分        *
%*                                                                 *
%* - Call procedures:                                              *
%*                          *
%*                                                                 *
%* - Called by :                                                   *
%*     SOLVE.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     Zixiong Lin, 2022.05.30                *
%*                                                                 *
%* *****************************************************************
function DynamicsSolve()

global cdata;
global sdata;

DSMETHOD = cdata.DSMETH; % 时间积分方法
switch DSMETHOD
    case 1
        FineInt(); % 精细积分法
    case 2
        ModalSup(); % 模态叠加法
    case 3
        Verlet(); % Verlet法
    case 4
        GenAlpha() % 广义alpha法
end

end