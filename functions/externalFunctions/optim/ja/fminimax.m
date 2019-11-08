% FMINIMAX   ���ϐ��֐��̃~�j�}�b�N�X���
%
% FMINIMAX �́A���̌`���̖��������܂��B
% {FUN(X)} �� X �Ɋւ��čŏ������܂��B�����ŁAFUN �� X �́A�x�N�g���A�܂���
% �s��ł��B
%
% X = FMINIMAX(FUN,X0) �́AX0 �������l�Ƃ��āA�֐� FUN �̃~�j�}�b�N�X�� 
% X �����߂܂��BFUN �́AX ����͂Ƃ��āAX �Ōv�Z�����x�N�g�� (�s��) ��
% �֐��l F ���o�͂��܂��BX0 �́A�X�J���A�x�N�g���A�܂��́A�s��ł��B
%
% X = FMINIMAX(FUN,X0,A,B) �́A���`�s�������� A*X <= B �̂��ƂŃ~�j�}�b�N�X
% ���������܂��B
%
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq) �́A���`�������� Aeq*X = Beq ���l������
% �~�j�}�b�N�X���������܂� (�s�������񂪂Ȃ��ꍇ�ɂ́AA=[],B=[] �Ɛݒ�
% ���܂�)�B
% 
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB) �́A�݌v�ϐ� X �̏㉺���͈̔͂�
% �^���邱�Ƃ��\�ł��B�����ŁALB <= X <= UB �͈͓��ōœK�����T������܂��B
% �͈͂̐��񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ肵�Ă��������BX(i) �ɉ�����
% �Ȃ��ꍇ�ALB(i) = -Inf �Ɛݒ肵�AX(i) �ɏ�����Ȃ��ꍇ�AUB(i) = Inf ��
% �ݒ肵�Ă��������B
% 
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) �́ANONLCON (�ʏ�� 
% M-�t�@�C�� NONLINCON.m) �Œ�`���ꂽ����̂��ƂŁA�~�j�}�b�N�X����
% �����܂��B�����Ŋ֐� NONLCON �́Afeval: [C, Ceq] = feval(NONLCON,X) ��
% �悤�ɌĂяo�����ꍇ�A���ꂼ��A����`�̕s��������Ɠ��������\�킷 
% C �� Ceq �x�N�g����Ԃ��܂��BFGOALATTAIN �́AC(X)< = 0 �� Ceq(X) = 0 ��
% �Ȃ�悤�ɍœK�����܂��B
%
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)  �́A�I�v�V���� 
% (OPTIONS) ��ݒ肵�Ď��s�ł��܂��B�����ŁAOPTIONS �� OPTIMSET �֐���
% �ݒ�ł���\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, 
% TolX, TolFun, TolCon, DerivativeCheck, FunValCheck, GradObj, GradConstr, 
% MaxFunEvals, MaxIter, MeritFunction, MinAbsMax, Diagnostics, DiffMinChange, 
% DiffMaxChange, OutputFcn, TypicalX �̃I�v�V�������g�p�ł��܂��BGradObj 
% �I�v�V�������g���āAFUN ���Ăяo���A2 �Ԗڂ̏o�͈��� G �ɓ_ X �ł�
% �Δ����W�� df/dX ��ݒ�ł��܂��B[F,G] = feval(FUN,X) �̌`���ŌĂяo���܂��B
% GradConstr �I�v�V�������g���āA���̂悤�ɁA4�̏o�͈����ŌĂяo����� 
% NONLCON ���w��ł��܂��B[C,Ceq,G,C,GCeq] = feval(NONLINCON,X) �̌`����
% �Ăяo���܂��B�����ŁAGC �́A�s��������x�N�g�� C �̕Δ����W���AGCeq �́A
% ��������x�N�g�� Ceq �̕Δ����W���ł��B�I�v�V������ݒ肵�Ȃ��ꍇ�́A
% OPTIONS = [] ���g�p���Ă��������B
%
% [X,FVAL] = FMINIMAX(FUN,X0,...) �́A�� X �ł̖ړI�֐� FVAL = feval(FUN,X) 
% �̒l���o�͂��܂��B
%
% [X,FVAL,MAXFVAL] = FMINIMAX(FUN,X0,...) �́A�� X �ł� MAXFVAL = max { FUN(X) } 
% ���o�͂��܂��B
%
% [X,FVAL,MAXFVAL,EXITFLAG] = FMINIMAX(FUN,X0,...) �́AFMINIMAX �̏I���󋵂�
% ���������� EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A
% ���̒ʂ�ł��B
%
%   1  1 ���̍œK�������w�肵�����e�͈͂𖞑����Ă��邱�Ƃ������܂��B
%   4  �T�������̑傫�����w�肵�����e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   5  �T�������ɉ������֐��̌��z�̑傫�������e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -1  �œK�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  ����������Ȃ��������Ƃ������܂��B
%   
% [X,FVAL,MAXFVAL,EXITFLAG,OUTPUT] = FMINIMAX(FUN,X0,...) �́A���s����
% �J��Ԃ��� OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A�ŏI
% �X�e�b�v�̃m���� OUTPUT.stepsize�A�ŏI�̃��C���T���̃X�e�b�v�� 
% OUTPUT.lssteplength�A�g�p�����A���S���Y�� OUTPUT.algorithm�A1 ���̍œK�� 
% OUTPUT.firstorderopt�A�I�����b�Z�[�W  OUTPUT.message �����\���� 
% OUTPUT ���o�͂��܂��B
%
% [X,FVAL,MAXFVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINIMAX(FUN,X0,...)  �́A
% ���ł̃��O�����W�F�搔 LAMBDA ���o�͂��܂��BLAMBDA �\���̂́A
% LAMBDA.lower �� LB ���ALAMBDA.upper �� UB ���ALAMBDA.ineqlin �ɐ��`
% �s�������ALAMBDA.eqlin �ɐ��`�������ALAMBDA.ineqnonlin �ɔ���`�s�������A
% LAMBDA.eqnonlin �ɔ���`������ݒ肵�܂��B
%
% ��
% FUN �́A@ ���g���Đݒ肷�邱�Ƃ��ł��܂��B
%
%        x = fminimax(@myfun,[2 3 4])
%
% �����ŁAmyfun �́A���̂悤�ɕ\����� MATLAB �֐��ł��B
% 
%   function F = myfun(x)
%   F = cos(x);
% 
% FUN �́A�����֐��Ƃ��Ă��\���ł��܂��B
% 
%       x = fminimax(@(x) sin(3*x),[2 5])
%
% FUN ���p�����[�^�����ꂽ�ꍇ�A���Ɉˑ������p�����[�^���w�肵�Ė���
% �֐����g�p�ł��܂��B2 �Ԗڂ̈��� c �Ńp�����[�^�����ꂽ�֐� myfun ��
% �ŏ�������Ɖ��肵�܂��B�����ŁAmfun �͂��̂悤�� M-�t�@�C���֐��ł��B
%
%       function F = myfun(x,c)
%       F = [x(1)^2 + c*x(2)^2;
%            x(2) - x(1)];
%
% c �̎w��l�ɑ΂��čœK�����邽�߂ɂ́A�ŏ��� c �ɐݒ肵�܂��B���ɁA
% x �������Ƃ��閳���֐� myfun ���`���A�ŏI�I�ɁAFMINIMAX �ɓn���܂��B
%
%       c = 2; % �ŏ��Ƀp�����[�^���`���܂��B
%       x = fminimax(@(x) myfun(x,c),[1;1])
%
%   �Q�l OPTIMSET, @, INLINE, FGOALATTAIN, LSQNONLIN.


%   Copyright 1990-2006 The MathWorks, Inc.
