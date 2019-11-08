% FMINCON  ����t�����ϐ��֐��̍ŏ���
%
% FMINCON �́A���̌`���̖��������܂��B
% ���̐���̂��ƂŁAF(X) ���AX �Ɋւ��čŏ������܂��B
% �@�@�@�@�@�@A*X  <= B, Aeq*X  = Beq (���`����)
%       �@�@�@C(X) <= 0, Ceq(X) = 0   (����`����)
%       �@�@�@LB <= X <= UB            
% 
% X = FMINCON(FUN,X0,A,B) �́A���� X0 �l�ŁA���`�s�������� A*X <= B ��
% ���ƂŁA�֐� FUN ���ŏ������� X �����߂܂��BFUN �́AX ����͂Ƃ��āA
% X �Ōv�Z�����X�J���̊֐��l���o�͂��܂��BX0 �́A�X�J���A�x�N�g���A
% �܂��́A�s��ł��B
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq) �́AA*X <= B �Ɠ��l�ɐ��`�������� 
% Aeq*X = Beq ���l������ FUN ���ŏ������܂� (�s�������񂪂Ȃ��ꍇ�́A
% A=[], B=[] �Ɛݒ肵�܂�)�B
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB) �́A�݌v�ϐ� X �̏㉺���͈̔͂�
% �^���邱�Ƃ��\�ł��B���̏ꍇ�ALB <= X <= UB �͈͓̔��ōœK�����T��
% ����܂��B�͈͂̐��񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ肵�Ă��������B
% X(i) �ɉ������Ȃ��ꍇ�ALB(i) = -Inf �Ɛݒ肵�AX(i) �ɏ�����Ȃ��ꍇ�A
% UB(i) = Inf �Ɛݒ肵�Ă��������B
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) �́ANONLCON �Œ�`���ꂽ
% ����̂��Ƃōŏ������܂��B�����Ŋ֐� NONLCON �́AX ����͂Ƃ��āA
% ����`�s�����Ɣ���`������\�킷 C �� Ceq �x�N�g����Ԃ��܂��BFMINCON �́A
% C(X) <= 0 �� Ceq(X) = 0 �ƂȂ�  FUN ���ŏ������܂��B(�͈͂̐��񂪂Ȃ��ꍇ�A
% LB=[] �� UB=[] �Ƃ��Đݒ肵�Ă��������B
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) �́A�I�v�V���� 
% (OPTIONS) ��ݒ肵�Ď��s�ł��܂��B�����ŁAOPTIONS �� OPTIMSET �֐���
% �ݒ�ł���\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������B
% Display, TolX, TolFun, TolCon, DerivativeCheck, Diagnostics, FunValCheck, 
% GradObj, GradConstr, Hessian, MaxFunEvals, MaxIter, DiffMinChange, 
% DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX, 
% Hessian, HessMult, HessPattern, OutputFcn �I�v�V�������g�p�ł��܂��BGradObj 
% �I�v�V�������g���āA2�Ԗڂ̏o�͈��� G �ɓ_ X �ł̊֐��̕Δ����W�� df/dX ��
% �o�͂��� FUN ��ݒ肵�܂��B�I�v�V���� Hessian ���g���āAFUN ���Ăяo���A
% 3�Ԗڂ̈��� H �ɁA�_ X �ł̃w�b�Z�s�� (2�K����) ��ݒ�ł��܂��BHessian �́A
% ��K�͖��݂̂Ɏg�p����A���C���T���@�ł͎g�p����܂���BGradConstr 
% �I�v�V�������g���āANONLCON ���A3 �Ԗڂ̏o�͈����� GC �A4 �Ԗڂ̈����� 
% GCeq �����悤�ɐݒ肵�܂��B�����ŁAGC �͕s��������x�N�g�� C �̕Δ���
% �W���AGCeq �́A��������x�N�g�� Ceq �̕Δ����W���ł��B�I�v�V������ݒ�
% ���Ȃ��ꍇ�́AOPTIONS = [] ���g�p���Ă��������B
% 
% [X,FVAL] = FMINCON(FUN,X0,...) �́A�� X �ł̖ړI�֐� FUN �̒l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG] = FMINCON(FUN,X0,...) �́AFMINCON �̏I���󋵂�����
% ������ EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A
% ���̒ʂ�ł��B
%
%   1  1 ���̍œK�������w�肵�����e�͈͂𖞑����Ă��邱�Ƃ������܂��B
%   2  X �̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   3  �ړI�֐��l�̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   4  �T�������̑傫�����w�肵�����e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   5  �T�������ɉ������֐��̌��z�̑傫�������e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ�
%      �����܂��B
%  -1  �œK�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  ����������Ȃ��������Ƃ������܂��B
% 
% [X,FVAL,EXITFLAG,OUTPUT] = FMINCON(FUN,X0,...) �́A
% �J��Ԃ��� OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A�ŏI
% �X�e�b�v�̃m���� OUTPUT.stepsize�A�g�p�����A���S���Y�� OUTPUT.algorithm�A
% 1 ���̍œK�� OUTPUT.firstorderopt�A�I�����b�Z�[�W OUTPUT.message ������
% �\���� OUTPUT ���o�͂��܂��B���K�̓A���S���Y���́A�ŏI�̃��C���T����
% �X�e�b�v�� OUTPUT.lssteplength ��Ԃ��A��K�̓A���S���Y���́A�������z
% �J��Ԃ��� OUTPUT.cgiterations ���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINCON(FUN,X0,...)  �́A���ł� 
% ���O�����W�F�搔 LAMBDA ���o�͂��܂��BLAMBDA �\���̂́ALAMBDA.lower �� 
% LB ���ALAMBDA.upper �� UB ���ALAMBDA.ineqlin �ɐ��`�s�������ALAMBDA.eqlin 
% �ɐ��`�������ALAMBDA.ineqnonlin �ɔ���`�s�������ALAMBDA.eqnonlin ��
% ����`������ݒ肵�܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD] = FMINCON(FUN,X0,...) �́A�� X ��
% �̊֐� FUN �̌��z�l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = FMINCON(FUN,X0,...) ��
% �� X �ł̊֐� FUN �̃w�b�Z�s����o�͂��܂��B
%�@
% ��F
% FUN �́A@ ���g���Đݒ肷�邱�Ƃ��ł��܂� : X = fmincon(@hump,...) 
% ���̏ꍇ�AF = humps(X) �́AX �ł� HUMPS �֐��̃X�J���֐��l F ���o�͂��܂��B
% 
% FUN �́A�����֐��Ƃ��Ă��\���ł��܂��B
% 
%        X = fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
% 
% �́AX = [0;0] ���o�͂��܂��B
% 
% FUN �܂��� NONLCON ���p�����[�^�����ꂽ�ꍇ�A���Ɉˑ������p�����[�^��
% �w�肵�Ė����֐����g�p�ł��܂��B�֐� myfun ���ŏ�������Ɖ��肵�A�܂�
% ����`������� mycon ���l�����܂��B�����ŁA������ 2 �̊֐��� 
% 2 �Ԗڂ̈��� a1 �� a2 �ł��ꂼ��p�����[�^������܂��B�����ŁAmyfun �� 
% mycon �͂��̂悤�� M-�t�@�C���֐��ł��B
%
%        function f = myfun(x,a1)
%        f = x(1)^2 + a1*x(2)^2;
%
% �����āA
%
%        function [c,ceq] = mycon(x,a2)
%        c = a2/x(1) - x(2);
%        ceq = [];
%
% �w�肳�ꂽ a1 �� a2 �ɑ΂��čœK�����s�Ȃ��ɂ́A�܂��A�����̃p�����[�^��
% �ݒ肵�܂��B�����āAa1 �� a2 �� 1 �̈����Ƃ��閳���֐��� 2 ��`���܂��B
% �ŏI�I�ɁA�����̖����֐��� FMINCON �ɓn���܂��B
%
%        a1 = 2; a2 = 1.5; % �ŏ��Ƀp�����[�^���`���܂��B
%        options = optimset('LargeScale','off'); % ���K�̓A���S���Y�������s�����܂��B
%        x = fmincon(@(x)myfun(x,a1),[1;2],[],[],[],[],[],[],@(x)mycon(x,a2),options)
%
%   �Q�l OPTIMSET, FMINUNC, FMINBND, FMINSEARCH, @, FUNCTION_HANDLE.


%   Copyright 1990-2006 The MathWorks, Inc.
