% FGOALATTAIN   ���ړI�S�[�����B�œK�����
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT) �́AX ��ω������邱�Ƃɂ��A�֐� 
% FUN �Őݒ肳���ړI�֐� (F) ���A�S�[���iGOAL�j�ɓ��B�����܂��B�S�[���́A
% WEIGHT �ɏ]���āA�d�ݕt�����܂��B������s�Ȃ��Ƃ��ɁA���̌`����
% ����`�v����������܂��B
% 
% LAMBDA :  F(X)-WEIGHT.*LAMBDA<=GOAL �� X, LAMBDA �Ɋւ��čŏ������܂��B
%           
% �֐� FUN �́AX ����͂Ƃ��āAX �Ōv�Z�����x�N�g�� (�s��) �̖ړI�l 
% F ���o�͂��܂��BX0 �́A�X�J���A�x�N�g���A�܂��͍s��ł��B
% 
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B) �́A���`�s�������� A*X <= B ��
% ���ƂŁA�S�[�����B���������܂��B
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq) �́A���`�������� 
% Aeq*X = Beq ���l�����āA�S�[�����B���������܂��B
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB)  �́A�݌v�ϐ� X 
% �̏㉺����^���邱�Ƃ��\�ł��B���̏ꍇ�́ALB <= X <= UB �͈͓̔���
% �œK�����T������܂��B�͈͂̐ݒ���s�Ȃ�Ȃ��ꍇ�́ALB �� UB �ɋ�s���
% �ݒ肵�Ă��������B�܂��A X(i) �ɉ������Ȃ��ꍇ�ALB(i) = -Inf �Ƃ��A
% X(i) �ɏ�����Ȃ��ꍇ�AUB(i) = Inf �Ɛݒ肵�Ă��������B
% 
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON) �́A
% NONLCON (�ʏ�́AM-�t�@�C�� NONLCON.m) �Œ�`��������̂��ƂŁA�S�[��
% ���B���������܂��B�����Ŋ֐� NONLCON �́A
% feval: [C, Ceq] = feval(NONLCON,X) �̂悤�ɌĂяo�����ꍇ�A���ꂼ��A
% ����`�̕s��������Ɠ��������\�킷 C �� Ceq �x�N�g����Ԃ��܂��B
% FGOALATTAIN �́AC(X)< = 0 �� Ceq(X) = 0 �ɂȂ�悤�ɍœK�����܂��B
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) 
% �́A�I�v�V���� (OPTIONS) ��ݒ肵�Ď��s�ł��܂��B�����ŁAOPTIONS �́A
% OPTIMSET �֐��Őݒ�ł���\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������B
% Display, TolX, TolFun, TolCon, DerivativeCheck, FunValCheck, GradObj, 
% GradConstr, MaxFunEvals, MaxIter, MeritFunction, GoalsExactAchieve, 
% Diagnostics, DiffMinChange, DiffMaxChange, OutputFcn, TypicalX ��
% �I�v�V�������g�p�ł��܂��BGradObj �I�v�V�������g���āAFUN ���Ăяo���A
% 2 �Ԗڂ̏o�͈��� G �� X �ł̕Δ����W�� df/dX ��ݒ�ł��܂��B
% [F,G] = feval(FUN,X) �̌`���ŌĂяo���܂��BGradConstr �I�v�V�������g���āA
% ���̂悤�ɁA4 �̏o�͈����ŌĂяo����� NONLCON ���w��ł��܂��B
% [C,Ceq,GC,GCeq] = feval(NONLCON,X) �̌`���ŌĂяo���܂��B�����ŁAGC �́A
% �s��������x�N�g�� C �̕Δ����W���AGCeq �́A��������x�N�g�� Ceq ��
% �Δ����W���ł��B�I�v�V������ݒ肵�Ȃ��ꍇ�́AOPTIONS = [] ���g�p����
% ���������B
%
% [X,FVAL] = FGOALATTAIN(FUN,X0,...) �́A�� X �ł̖ړI�֐� FUN �̒l��
% �o�͂��܂��B
%
% [X,FVAL,ATTAINFACTOR] = FGOALATTAIN(FUN,X0,...) �́A�� X �ł̓��B�t�@�N�^
% ���o�͂��܂��BATTAINFACTOR �����̏ꍇ�A�S�[���͉ߓ��B�ɂȂ�܂��B�܂��A
% ���̏ꍇ�A�����B�ɂȂ�܂��B
%
% [X,FVAL,ATTAINFACTOR,EXITFLAG] = FGOALATTAIN(FUN,X0,...) �́AFGOALATTAIN ��
% �I���󋵂����������� EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����
% �I���󋵂́A���̒ʂ�ł��B
%
%   1  FGOALATTAIN �́A�� X �Ɏ������Ă��邱�Ƃ������܂��B
%   4  �T�������̑傫�����w�肵�����e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   5  �T�������ɉ������֐��̌��z�̑傫�������e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -1  �œK�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  ����������Ȃ��������Ƃ������܂��B
%   
% [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT] = FGOALATTAIN(FUN,X0,...) �́A
% �J��Ԃ��� OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A
% �ŏI�X�e�b�v�̃m���� OUTPUT.stepsize�A�ŏI�̃��C���T���̃X�e�b�v�� 
% OUTPUT.lssteplength�A�g�p�����A���S���Y�� OUTPUT.algorithm�A1 ���̍œK�� 
% OUTPUT.firstorderopt�A�I�����b�Z�[�W  OUTPUT.message �����\���� 
% OUTPUT ���o�͂��܂��B
% 
% [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT,LAMBDA] = FGOALATTAIN(FUN,X0,...)
% �́A���ł̃��O�����W�F�搔 LAMBDA ���o�͂��܂��BLAMBDA �\���̂́A
% LAMBDA.lower �� LB ���ALAMBDA.upper �� UB ���ALAMBDA.ineqlin �ɐ��`
% �s�������ALAMBDA.eqlin �ɐ��`�������ALAMBDA.ineqnonlin �ɔ���`�s�������A
% LAMBDA.eqnonlin �ɔ���`������ݒ肵�܂��B
%
% �ڍׂ́AM-�t�@�C�� FGOALATTAIN.M ���Q�Ƃ��Ă��������B
%
% �Q�l OPTIMSET, OPTIMGET.


%   Copyright 1990-2006 The MathWorks, Inc.
