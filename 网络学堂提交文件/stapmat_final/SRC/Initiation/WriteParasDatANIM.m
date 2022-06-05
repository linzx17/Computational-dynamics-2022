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
%*     Engineering, Tsinghua University, 2021.05.12                *
%*                                                                 *
%* *****************************************************************

function WriteParasDatANIM()
global cdata;

% Open file
IDAT_ANIM = cdata.IDAT_ANIM;
fprintf(IDAT_ANIM,'TITLE = "Example : %s"\n',cdata.HED);
fprintf(IDAT_ANIM,'VARIABLES= "X"   "Y"   "Z" \n\n');

end

