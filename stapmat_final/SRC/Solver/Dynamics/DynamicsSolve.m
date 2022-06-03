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
        PreciseIntegral(); % 精细积分法最终版
    case 2
        ModalSup(); % 模态叠加法
    case 3
        Verlet(); % 速度Verlet法
    case 4
        GenAlpha(); % 广义alpha法
    case 5
        ModalSup(); % 模态叠加法
        Verlet(); % 速度Verlet法
        GenAlpha() % 广义alpha法
        PreciseIntegral_quick(); % 惊喜积分法快速版
        PreciseIntegral(); % 精细积分法最终版
    case 6
        PreciseInt(); % 精细积分法（第一版精细积分）
    case 7
        PreciseIntSPARSE(); % 具有矩阵稀疏化的精细积分法（第二版精细积分）
    case 8
        PreciseIntegral_quick(); % 精细积分法快速版
    otherwise
        fprintf('不进行时间积分计算');
end

end