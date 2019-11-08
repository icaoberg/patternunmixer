% LSQNONLIN   ����`�ŏ������
%
% LSQNONLIN �́A���̌`���̖��������܂��B
% sum {FUN(X).^2} ���AX �Ɋւ��čŏ������܂��B�����ŁAX �Ɗ֐� FUN ��
% �Ԃ����l�́A�x�N�g���A�܂��́A�s��ł��B
%
% X = LSQNONLIN(FUN,X0) �́AX0 �������l�Ƃ��āA�֐� FUN �ɋL�q�����֐�
% �̓��a���ŏ������� X �����߂܂��B�֐� FUN �́AX ����͂Ƃ��āAX ��
% �v�Z�����x�N�g�� (�܂��͍s��) �̊֐��l F ���o�͂��܂��B���� �F FUN �́A
% FUN(X) ��Ԃ��A���a sum(FUN(X).^2) �͕Ԃ��܂��� (���̏����́A�A���S���Y��
% ���ŉA�I�ɍs�Ȃ��܂�)�B
%
% X = LSQNONLIN(FUN,X0,LB,UB) �́A�݌v�ϐ� X �̏㉺���͈̔͂�^���邱�Ƃ�
% �\�ł��B���̏ꍇ�ALB <= X <= UB �͈̔͂ōœK�����T������܂��B�͈͂�
% ���񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ肵�Ă��������BX(i) �ɉ������Ȃ��ꍇ�A
% LB(i) = -Inf �Ɛݒ肵�AX(i) �ɏ�����Ȃ��ꍇ�AUB(i) = Inf �Ɛݒ肵�Ă��������B
% 
% X = LSQNONLIN(FUN,X0,LB,UB,OPTIONS) �́A�I�v�V���� (OPTIONS) ��ݒ肵��
% ���s�ł��܂��B������ OPTIONS �� OPTIMSET �֐��Őݒ�ł���\���̂ł��B
% �ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, TolX, TolFun, DerivativeCheck, 
% Diagnostics, FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType, 
% LevenbergMarquardt, MaxFunEvals, MaxIter, DiffMinChange, DiffMaxChange, 
% LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX, OutputFcn 
% �̃I�v�V�������g�p�ł��܂��B�I�v�V���� Jacobian ���g���āAFUN ���Ăяo���A
% 2 �Ԗڂ̏o�͈��� J �ɁA�_ X �ł̃��R�r�s���ݒ�ł��܂��BX ������ n ��
% �Ƃ��ɁAFUN �� m �v�f�̃x�N�g�� F ��Ԃ��ꍇ�AJ �́Am �s n ��̍s���
% �Ȃ�܂��B�����ŁAJ(i,j) �́AF(i) �� x(j) �ɂ��Δ����W���ł� (���R�r�s�� 
% J �́AF �̌��z��]�u�������̂ł�)�B
%
% [X,RESNORM] = LSQNONLIN(FUN,X0,...) �́AX �ł̎c���� 2 �m���� sum(FUN(X).^2) 
% ���o�͂��܂��B
%
% [X,RESNORM,RESIDUAL] = LSQNONLIN(FUN,X0,...) �́A�� X �ł̎c���l 
% RESIDUAL = FUN(X) ��Ԃ��܂��B
%
% [X,RESNORM,RESIDUAL,EXITFLAG] = LSQNONLIN(FUN,X0,...) �́ALSQNONLIN ��
% �I���󋵂����������� EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����
% �I���󋵂́A���̒ʂ�ł��B
%
%   1  LSQNONLIN �́A�� X �Ɏ����������Ƃ������܂��B
%   2  X �̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   3  �ړI�֐��l�̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   4  �T�������̑傫�����w�肵�����e�덷��菬�������Ƃ������܂��B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -1  �A���S���Y�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  ���E���������Ă��邱�Ƃ������܂��B
%  -4  ���C���T���́A���݂̒T�������ɉ����ď\���Ɏc���������ł��Ȃ����Ƃ������܂��B
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQNONLIN(FUN,X0,...) �́A�J��
% �Ԃ��� OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A�g�p����
% �A���S���Y�� OUTPUT.algorithm�A(�g�p�����ꍇ) �������z�J��Ԃ��� 
% OUTPUT.cgiterations�A(�g�p�����ꍇ) 1 ���̍œK�� OUTPUT.firstorderopt�A
% �I�����b�Z�[�W OUTPUT.message �����\���� OUTPUT ���o�͂��܂��B
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQNONLIN(FUN,X0,...) �́A
% ���ł̃��O�����W�F�搔 LAMBDA ���o�͂��܂��BLAMBDA.lower �� LB ���A
% LAMBDA.upper ��UB ��ݒ肵�܂��B
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = LSQNONLIN(FUN,X0,...)
% �́AX �ł̊֐� FUN �̃��R�r�s����o�͂��܂��B
%
% ��
% FUN �́A@ ���g���Đݒ肷�邱�Ƃ��ł��܂��B
%        x = lsqnonlin(@myfun,[2 3 4])
%
% �����ŁAmyfun �́A���̂悤�ɕ\����� MATLAB �֐��ł��B
%
%       function F = myfun(x)
%       F = sin(x);
%
% FUN �́A�����֐��Ƃ��Ă��\���ł��܂��B
%
%       x = lsqnonlin(@(x) sin(3*x),[1 4])
%
% FUN �܂��� NONLCON ���p�����[�^�����ꂽ�ꍇ���Ɉˑ������p�����[�^��
% �w�肵�Ė����֐����g�p�ł��܂��B2 �Ԗڂ̈��� c �Ńp�����[�^�����ꂽ�֐� 
% myfun ���ŏ�������Ɖ��肵�܂��B�����ŁAmfun �͂��̂悤�� M-�t�@�C��
% �֐��ł��B
%
%       function F = myfun(x,c)
%       F = [ 2*x(1) - exp(c*x(1))
%             -x(1) - exp(c*x(2))
%             x(1) - x(2) ];
%
% c �̎w��l�ɑ΂��čœK�����邽�߂ɂ́A�ŏ��� c �ɐݒ肵�܂��B���ɁA
% x �������Ƃ��閳���֐� myfun ���`���A�ŏI�I�ɁALSQNONLIN �ɓn���܂��B
%
%       c = -1; % �ŏ��Ƀp�����[�^���`���܂��B
%       x = lsqnonlin(@(x) myfun(x,c),[1;1])
%
%   �Q�l OPTIMSET, LSQCURVEFIT, FSOLVE, @, INLINE.


%   Copyright 1990-2006 The MathWorks, Inc.
