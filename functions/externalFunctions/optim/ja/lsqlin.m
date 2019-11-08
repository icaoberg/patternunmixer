% LSQLIN   ����t�����`�ŏ������
%
% X = LSQLIN(C,d,A,b) �́A���̌`���̕s�����̐�������t���ŏ�������
% �����܂��B
%
% A*x <=  b �̐���̂��ƂŁA0.5*(NORM(C*x-d)).^2 �� x �Ɋւ��čŏ������܂��B
% �����ŁAC �́Am �s n ��̍s��ł��B
%
% X = LSQLIN(C,d,A,b,Aeq,beq) �́A���̌`���� (�����̐���t��) �ŏ����
% ���������܂��B
%
% A*x <=  b �� Aeq*X = Beq ���l�����āA0.5*(NORM(C*x-d)).^2 �� x �Ɋւ��āA
% �ŏ������܂��B
%
% X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB) �́A�݌v�ϐ� X �̏㉺���͈̔͂�^����
% ���Ƃ��\�ł��B���̏ꍇ�ALB <= X <= UB �͈͓̔��ōœK�����T������܂��B
% �͈͂̐��񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ肵�Ă��������BX(i) ��
% �������Ȃ��ꍇ�ALB(i) = -Inf �Ɛݒ肵�AX(i) �ɏ�����Ȃ��ꍇ�AUB(i) = Inf 
% �Ɛݒ肵�Ă��������B
%
% X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0) �́AX0 �������l�Ƃ��܂��B
%
% X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0,OPTIONS) �́A�I�v�V���� (OPTIONS) 
% ��ݒ肵�Ď��s�ł��܂��B������ OPTIONS �� OPTIMSET �֐��Őݒ�ł���
% �\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, Diagnostics, 
% TolFun, LargeScale, MaxIter, JacobMult, PrecondBandWidth, TypicalX, 
% TolPCG, MaxPCGIter �̃I�v�V�������g�p�ł��܂��B'final' �� 'off' �݂̂��A
% Display �ɑ΂��Ďg���܂� ('iter' �͎g�p�ł��܂���)�B
%
% X=LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0,OPTIONS,P1,P2,...) �́A
% OPTIMSET('JacobMult',JMFUN) ���ݒ肳�ꂽ�ꍇ�A���Ɉˑ�����p�����[�^ 
% P1,P2,... �𒼐� JMFUN �֐��ɓn���܂��BJMFUN �́A���[�U�ɂ���Ďw��
% ���܂��B�f�t�H���g�l���g�p����ɂ́AA, b, Aeq, beq, LB, UB, XO, 
% OPTIONS �ɋ�s���n���Ă��������B
%
% [X,RESNORM] = LSQLIN(C,d,A,b) �́A�c���̓�� 2 �m�����l norm(C*X-d)^2 ��
% �o�͂��܂��B
%
% [X,RESNORM,RESIDUAL] = LSQLIN(C,d,A,b) �́A�c�� C*X-d ���o�͂��܂��B
%
% [X,RESNORM,RESIDUAL,EXITFLAG] = LSQLIN(C,d,A,b) �́ALSQLIN �̏I���󋵂�
% ���������� EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A
% ���̒ʂ�ł��B
%
%   1  LSQLIN �́A�� X �Ɏ����������Ƃ������܂��B
%   3  �c���̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   0  �J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -2  ��肪�s���ł��邱�Ƃ������܂��B
%  -4  �������������߁A�œK���𒆒f�������Ƃ������܂��B
%  -7  �T�������̑傫�����������Ȃ肷���Ă��邱�Ƃ������܂��B����ȏ�
%      �v�Z���s�Ȃ����Ƃ��ł��܂���B
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQLIN(C,d,A,b) �́A�J��Ԃ� 
% OUTPUT.iterations�A�g�p�����A���S���Y�� OUTPUT.algorithm�A(�g�p�����ꍇ) 
% �������z�J��Ԃ��� OUTPUT.cgiterations�A(��K�͖��̂�) 1 ���̍œK�� 
% OUTPUT.firstorderopt�A�I�����b�Z�[�W OUTPUT.message �����\���� OUTPUT ��
% �o�͂��܂��B
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQLIN(C,d,A,b)  �́A
% �œK���ł̃��O�����W�F�搔 LAMBDA ���o�͂��܂��BLAMBDA.ineqyalities �ɐ��`
% �s���� C ���ALAMBDA.eqlin �ɐ��`���� Ceq ���ALAMBDA.lower �� LB ���A
% LAMBDA.upper �� UB ��ݒ肵�܂��B


%   Copyright 1990-2006 The MathWorks, Inc.
