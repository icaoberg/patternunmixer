% FSOLVE   �ŏ����@���g�������ϐ�����`�������̋���
%
% FSOLVE �́A���̌`���̖��������܂��B
% 
% F(X) = 0    �����ŁAF �� X �́A�x�N�g���܂��͍s��ł��B  
%
% X = FSOLVE(FUN,X0) �́A�����l���s�� X0 �Ƃ��āA�֐� FUN �ɂ��\�����
% �������������܂��B�֐� FUN �́AX ����͂Ƃ��āAX �Ōv�Z�����x�N�g�� 
% (�s��) �̊֐��l F ���o�͂��܂��B
%
% X = FSOLVE(FUN,X0,OPTIONS) �́A�I�v�V���� (OPTIONS) ��ݒ肵�Ď��s
% �ł��܂��B������ OPTIONS �� OPTIMSET �֐��Őݒ�ł���\���̂ł��B�ڍׂ́A
% OPTIMSET ���Q�Ƃ��Ă��������BDisplay, TolX, TolFun, DerivativeCheck, 
% Diagnostics, FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType,
% NonlEqnAlgorithm, MaxFunEvals, MaxIter, OutputFcn, DiffMinChange, DiffMaxChange, 
% LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX �̃I�v�V������
% �g�p�ł��܂��B�I�v�V���� Jacobian ���g���āAFUN ���Ăяo���A2 �Ԗڂ�
% �o�͈��� J �ɁA�_ X �ł̃��R�r�s���ݒ�ł��܂��BFUN �́AX ������ 
% n �̂Ƃ��ɁAm �v�f�̃x�N�g�� F ���o�͂���ꍇ�AJ �́Am �s n ��̍s���
% �Ȃ�܂��B�����ŁAJ(i,j) �́AF(i) �� x(j) �ɂ��Δ����W���ł� 
% (���R�r�s�� J �́AF �̌��z��]�u�������̂ł�)�B
%
% [X,FVAL] = FSOLVE(FUN,X0,...) �́AX �ł̕����� FUN �̒l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG] = FSOLVE(FUN,X0,...) �́AFSOLVE �̏I���󋵂�����
% ������ EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A
% ���̒ʂ�ł��B
%
%   1  FSOLVE �́A�� X �Ɏ����������Ƃ������܂��B
%   2  X �̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   3  �ړI�֐��l�̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   4  �T�������̑傫�����w�肵�����e�덷��菬�������Ƃ������܂��B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -1  �œK�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  �A���S���Y�������łȂ��_�Ɏ������Ă��邱�Ƃ������܂��B
%  -3  �M���̈攼�a���������Ȃ肷���Ă��邱�Ƃ������܂��B
%  -4  ���C���T���́A���݂̒T�������ɉ����ď\���Ɏc���������ł��Ȃ����Ƃ������܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT] = FSOLVE(FUN,X0,...) �́A�J��Ԃ��� 
% OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A�g�p�����A���S���Y�� 
% OUTPUT.algorithm�A(�g�p�����ꍇ) �������z�J��Ԃ��� OUTPUT.cgiterations�A
% (�g�p�����ꍇ) 1 ���̍œK�� OUTPUT.firstorderopt�A�I�����b�Z�[�W 
% OUTPUT.message �����\���� OUTPUT ���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,JACOB] = FSOLVE(FUN,X0,...)  �́AX �ł̊֐�
% FUN �̃��R�r�s����o�͂��܂��B
% 
% ��
% FUN �́A@ ���g���Đݒ肷�邱�Ƃ��ł��܂��B
%
%        x = fsolve(@myfun,[2 3 4],optimset('Display','iter'))
%
% �����ŁAmyfun �́A���̂悤�ɕ\����� MATLAB �֐��ł��B
%
%       function F = myfun(x)
%       F = sin(x);
%
% FUN �́A�����֐��Ƃ��Ă��\���ł��܂��B
%
%       x = fsolve(@(x) sin(3*x),[1 4],optimset('Display','off'))
%
% FUN ���p�����[�^�����ꂽ�ꍇ�A���Ɉˑ������p�����[�^���w�肵�Ė���
% �֐����g�p�ł��܂��B2 �Ԗڂ̈��� c �Ńp�����[�^�����ꂽ�֐� myfun ��
% �ŏ�������Ɖ��肵�܂��B�����ŁAmyfun �͂��̂悤��M-�t�@�C���֐��ł��B
%     
%       function F = myfun(x,c)
%       F = [ 2*x(1) - x(2) - exp(c*x(1))
%             -x(1) + 2*x(2) - exp(c*x(2))];
%
% c �̎w��l�ɑ΂��čœK�����邽�߂ɂ́A�ŏ��� c �ɐݒ肵�܂��B���ɁA
% x �������Ƃ��閳���֐� myfun ���`���A�ŏI�I�ɁAFSOLVE �ɓn���܂��B
%
%       c = -1; % �ŏ��Ƀp�����[�^���`���܂��B
%       x = fsolve(@(x) myfun(x,c),[-5;-5])
%
%   �Q�l OPTIMSET, LSQNONLIN, @, INLINE.


%   Copyright 1990-2006 The MathWorks, Inc.
