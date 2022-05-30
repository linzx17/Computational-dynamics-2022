% % 根据不同控制变量，调用不同的求解器进行求解
function SOLVE()

global cdata;
global sdata;

MODEX = cdata.MODEX;

if (1 == MODEX) || (2 == MODEX) || (3 == MODEX) %% 静力学求解
    StaticsSolve();
    if 3 == MODEX
        Eigenvalue(); %% 模态求解
    end
elseif (4 == MODEX)
    Eigenvalue(); %% 模态求解
    DynamicsSolve(); %% 动力学求解
end

end
