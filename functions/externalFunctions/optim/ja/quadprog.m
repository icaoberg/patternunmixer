% QUADPROG  �񎟌v��@
%
% X = QUADPROG(H,f,A,b) �́A���̌`���̓񎟌v����������܂��B
%
% A*x <=  b �̂��ƂŁA0.5*x'*H*x + f'*x �� x �Ɋւ��čŏ������܂��B
% 
% X = QUADPROG(H,f,A,b,Aeq,beq) �́A���`�������� Aeq*x = beq �̂��ƂŁA
% ��L�̖��������܂��B
%
% X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB) �́A�݌v�ϐ� X �̏㉺���͈̔͂�
% �^���邱�Ƃ��\�ł��B���̏ꍇ�ALB <= X <= UB �͈͓̔��ōœK����
% �T������܂��B�͈͂̐��񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ肵�Ă��������B
% X(i) �ɉ������Ȃ��ꍇ�ALB(i) = -Inf �Ɛݒ肵�AX(i) �ɏ�����Ȃ��ꍇ�A
% UB(i) = Inf �Ɛݒ肵�Ă��������B
%
% X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0) �́AX0 �������l�Ƃ��܂��B
%
% X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS)  �́A�I�v�V���� (OPTIONS) 
% ��ݒ肵�Ď��s�ł��܂��B������ OPTIONS �� OPTIMSET �֐��Őݒ�ł���\����
% �ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, Diagnostics, TolX, 
% TolFun, HessMult, LargeScale, MaxIter, PrecondBandWidth, TypicalX, 
% TolPCG, MaxPCGIter �̃I�v�V�������g�p�ł��܂��B'final' �� 'off' �݂̂�
% Display �ɑ΂��Ďg���܂� ('iter' �͎g�p�ł��܂���)�B
%
% X=QUADPROG(Hinfo,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS,P1,P2,...) �́A
% OPTIMSET('HessMult',HMFUN) ���ݒ肳�ꂽ�ꍇ�A���Ɉˑ�����p�����[�^ 
% P1, P2,... �𒼐� HMFUN �֐��ɓn���܂��BHMFUN �́A���[�U�ɂ���ė^��
% ���܂��B�f�t�H���g�l���g�p����ꍇ�́AA, b, Aeq, beq, LB, UB, XO, 
% OPTIONS �ɋ�s���n���Ă��������B
%
% [X,FVAL] = QUADPROG(H,f,A,b) �́AX �ł̖ړI�֐� FVAL = 0.5*X'*H*X + f'*X 
% �̒l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG] = QUADPROG(H,f,A,b) �́AQUADPROG �̏I���󋵂�����
% ������ EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A
% ���̒ʂ�ł��B
%
%   1  QUADPROG �́A�� X �Ɏ����������Ƃ������܂��B
%   3  �ړI�֐��l�̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   4  �Ǐ��ŏ��l�������������Ƃ������܂��B
%   0  �J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -2  ����������Ȃ��������Ƃ������܂��B
%  -3  ���ɐ��񂪂Ȃ����Ƃ������܂��B
%  -4  �������������߁A�œK���𒆒f�������Ƃ������܂��B
%  -7  �T�������̑傫�����������Ȃ肷���Ă��邱�Ƃ������܂��B����ȏ�
%      �v�Z���s�Ȃ����Ƃ��ł��܂���B
%
% [X,FVAL,EXITFLAG,OUTPUT] = QUADPROG(H,f,A,b) �́A�J��Ԃ��� 
% OUTPUT.iterations�A�g�p�����A���S���Y�� OUTPUT.algorithm�A(�g�p����
% �ꍇ) �������z�J��Ԃ��� OUTPUT.cgiterations �A(��K�͖@�̂�) 1 ����
% �œK�l OUTPUT.firstorderopt�A�I�����b�Z�[�W������ꍇ�� OUTPUT.message �� 
% ���\���� OUTPUT ���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = QUADPROG(H,f,A,b) �́A�œK���ł�
% ���O�����W�F�搔 LAMBDA ���o�͂��܂��B LAMBDA.ineqlin ��
% ���`�s���� A ���ALAMBDA.eqlin �ɐ��`���� Aeq ���ALAMBDA.lower �� 
% LB ���ALAMBDA.upper �� UB ��ݒ肵�Ă��܂��B


%   Copyright 1990-2006 The MathWorks, Inc. 
