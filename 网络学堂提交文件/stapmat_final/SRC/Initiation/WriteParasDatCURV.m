%* *****************************************************************
%* - Function of STAPMAT in initialization phase                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Write parameters to output file for postdeal for Tecplot                 *
%*                                                                 *
%* - Call procedures: None                                         *
%*                                                                 *
%* - Called by :                                                   *
%*     stapmat.m                                                   *
%*                                                                 *
%* - Programmed by:                                                *
%*     Wang Wanting                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2021.05.15                *
%*                                                                 *
%* *****************************************************************

function WriteParasDatCURV()
global cdata;

% Open file
IDAT_CURV = cdata.IDAT_CURV;
fprintf(IDAT_CURV,'TITLE = "Example : %s"\n',cdata.HED);
fprintf(IDAT_CURV,['VARIABLES="NODE_ID"       "U1"           "U2"           "U3"' ...
    '          "E11"          "E22"          "E33"          "E12"          "E13"          "E23"' ...
    '          "S11"          "S22"          "S33"          "S12"          "S13"          "S23"\n']);
end

