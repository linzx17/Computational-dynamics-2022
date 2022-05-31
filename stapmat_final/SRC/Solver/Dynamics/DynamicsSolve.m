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
        PreciseInt(); % 精细积分法
    case 2
        PreciseIntSPARSE(); % 具有矩阵稀疏化的精细积分法
    case 3
        ModalSup(); % 模态叠加法
    case 4
        Verlet(); % 速度Verlet法
    case 5
        GenAlpha() % 广义alpha法
    case 6
        PreciseIntgral(); % 精细积分法最终版
    case 7
        ModalSup(); % 模态叠加法
        Verlet(); % 速度Verlet法
        GenAlpha() % 广义alpha法
end

end