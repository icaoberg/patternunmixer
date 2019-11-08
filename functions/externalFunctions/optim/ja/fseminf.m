% FSEMINF   ����������t���œK�����
%
% FSEMINF �́A���̌`���̖��������܂��B
%
%         { F(x) | C(x) <= 0 , Ceq(X) = 0 , PHI(x,w) <= 0 } 
% ��ԓ��ɑ��݂��邷�ׂĂ� w �ɑ΂��āAF(x) �� x �Ɋւ��čŏ������܂��B
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON) �́AX0 �������l�Ƃ��āA�֐� SEMINFCON 
% �ŋK�肳�ꂽ NTHETA �̔���������̂��ƂŁA�֐� FUN �ŏ������� X ��
% ���߂܂��B�֐� FUN �́AX ���x�N�g�����͂Ƃ��āAX �Ōv�Z�����X�J����
% �֐��l F ���o�͂��܂��B�֐� SEMINFCON �́AX, S �����`���͂Ƃ��āA
% ����`�s��������x�N�g�� C�A����`��������x�N�g�� Ceq�A����сA����
% ��ԑS�̂ŕ]������� NTHETA �̔������s��������s�� PHI_1, PHI_2, ..., 
% PHI_NTHETA ���o�͂��܂��BS �́A��������T���v���Ԋu�ŁA�g�p���Ă����Ȃ��Ă�
% �\���܂���B
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B) �́A���`�s�������� A*X <=  B ��
% ���������悤�Ƃ��܂��B
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq) �́A���`�������� 
% Aeq*X = Beq �����������悤�Ƃ��܂� (�s���������݂��Ȃ��ꍇ�AA = [] �� 
% B = [] ��ݒ肵�Ă�������)�B
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq,LB,UB) �́A�݌v�ϐ� X 
% �̏㉺���͈̔͂�^���邱�Ƃ��\�ł��B���̏ꍇ�ALB <= X <= UB �͈͓̔���
% �œK�����T������܂��B�͈͂̐��񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ�
% ���Ă��������BX(i) �ɉ������Ȃ��ꍇ�ALB(i) = -Inf �Ɛݒ肵�AX(i) ��
% ������Ȃ��ꍇ�AUB(i) = Inf �Ɛݒ肵�Ă��������B
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq,LB,UB,OPTIONS) �́A
% �I�v�V���� (OPTIONS) ��ݒ肵�Ď��s�ł��܂��B�����ŁAOPTIONS �� OPTIMSET 
% �֐��Őݒ�ł���\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, 
% TolX, TolFun, TolCon, DerivativeCheck, Diagnostics, FunValCheck, GradObj, 
% MaxFunEvals, MaxIter, DiffMinChange, DiffMaxChange, OutputFcn, TypicalX 
% �̃I�v�V�������g�p�ł��܂��B�I�v�V���� GradObj ���g���āAFUN ���Ăяo���A
% 2 �Ԗڂ̏o�͈��� G �ɓ_ X �ł̕Δ����W�� df/dX ��ݒ�ł��܂��B
% [F,G] = feval(FUN,X) �̌`���ŌĂяo���܂��B
%
% [X,FVAL] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) �́A�� X �ł̖ړI�֐�
% FUN �̒l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) �́AFSEMINF ��
% �I���󋵂����������� EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����
% �I���󋵂́A���̒ʂ�ł��B
% 
%   1  FSEMINF �͉� X �Ɏ����������Ƃ������܂��B
%   4  �T�������̑傫�����w�肵�����e�덷��菬�����A����ᔽ�� options.TolCon 
%      ��菬�������Ƃ������܂��B
%   5  �T�������ɉ������֐��̌��z�̑傫�������e�͈͂�菬�����A����ᔽ�� 
%      options.TolCon ��菬�������Ƃ������܂��B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -1  �œK�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  ����������Ȃ��������Ƃ������܂��B
% 
% [X,FVAL,EXITFLAG,OUTPUT] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) �́A
% �J��Ԃ��� OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A
% �ŏI�X�e�b�v�̃m���� OUTPUT.stepsize�A�ŏI�̃��C���T���̃X�e�b�v�� 
% OUTPUT.lssteplength�A�g�p�����A���S���Y�� OUTPUT.algorithm�A1 ���̍œK�� 
% OUTPUT.firstorderopt�A�I�����b�Z�[�W  OUTPUT.message �����\���� 
% OUTPUT ���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...)
% �́A���ł̃��O�����W�F�搔���o�͂��܂��BLAMBDA �\���̂́ALAMBDA.lower
% �� LB ���ALAMBDA.upper �� UB ���ALAMBDA.ineqlin �ɐ��`�s�������ALAMBDA.eqlin 
% �ɐ��`�������ALAMBDA.ineqnonlin �ɔ���`�s�������ALAMBDA.eqnonlin ��
% ����`������ݒ肵�܂��B
% 
% ��
% FUN �� SEMINFCOM �́A@ ���g���āA�ݒ肷�邱�Ƃ��ł��܂��B
%
%        x = fseminf(@myfun,[2 3 4],3,@myseminfcon)
%
% �����ŁAmyfun �́A���̂悤�ɕ\�킳��� MATLAB �֐��ł��B
% 
%    function F = myfun(x)
%    F = x(1)*cos(x(2))+x(3)^3:
% 
% �܂��Amyseminfcon �́A���̂悤�ɕ\�킳��� MATLAB �֐��ł��B
% 
%       function [C,Ceq,PHI1,PHI2,PHI3,S] = myseminfcon(X,S)
%       C = [];     % C �� Ceq ���v�Z����R�[�h:
%                   % ���񂪂Ȃ��ꍇ�͋�s��ɂ��܂��B
%       Ceq = [];
%       if isnan(S(1,1))
%          S = [...] ; % S �� ntheta�s�~2��
%       end
%       PHI1 = ... ;       % PHI ���v�Z����R�[�h
%       PHI2 = ... ;
%       PHI3 = ... ;
% �@�@�@�@�@
%   �Q�l OPTIMSET, @, FGOALATTAIN, LSQNONLIN.


%   Copyright 1990-2006 The MathWorks, Inc.
