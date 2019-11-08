% FMINUNC   ���ϐ��֐��̍ŏ���
%
% X = FMINUNC(FUN,X0) �́AX0 �������l�Ƃ��āA�֐� FUN �̋Ǐ��ŏ��� X ��
% �����܂��BFUN �́AX ����͂Ƃ��āAX �Ōv�Z�����X�J���̊֐��l F ��
% �o�͂��܂��BX0 �́A�X�J���A�x�N�g���A�܂��́A�s��ł��B
% 
% X = FMINUNC(FUN,X0,OPTIONS) �́A�I�v�V���� (OPTIONS) ��ݒ肵�Ď��s
% �ł��܂��B������ OPTIONS �� OPTIMSET �֐��Őݒ�ł���\���̂ł��B
% �ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, TolX, TolFun, DerivativeCheck, 
% Diagnostics, FunValCheck, GradObj, HessPattern, Hessian, HessMult, HessUpdate, 
% InitialHessType, InitialHessMatrix, MaxFunEvals, MaxIter, DiffMinChange, 
% DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, OutputFcn, 
% TypicalX �̃I�v�V�������g�p�ł��܂��BGradObj �I�v�V�������g���āAFUN ��
% �Ăяo���A2 �Ԗڂ̏o�͈��� G �ɓ_ X �ł̕Δ����W�� df/dX ��ݒ�ł��܂��B
% �I�v�V���� Hessian ���g���āAFUN ���Ăяo���A3 �Ԗڂ̈��� H �ɓ_ X �ł�
% �w�b�Z�s�� (2 �K����) ��ݒ�ł��܂��BHessian �́A��K�͖��ł̂ݎg���A
% ���C���T���@�ł͎g���܂���B
%
% [X,FVAL] = FMINUNC(FUN,X0,...) �́A�� X �ł̖ړI�֐� FUN �̒l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG] = FMINUNC(FUN,X0,...) �́AFMINUNC �̏I���󋵂�����
% ������  EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A
% ���̒ʂ�ł��B
%
%   1  ���z�̑傫�����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   2  X �̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂��B
%   3  �ړI�֐��l�̕ω����w�肵�����e�͈͂�菬�������Ƃ������܂�
%      (��K�͖��̂�)�B
%   0  �֐��v�Z�̉񐔁A���邢�͌J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ�
%      �����܂��B
%  -1  �A���S���Y�����A�o�͊֐��ŏI�����Ă��邱�Ƃ������܂��B
%  -2  ���C���T�����A���݂̒T�������ɉ����Ď󂯓���\�ȓ_������
%      ���Ȃ����Ƃ������܂� (���K�͖��̂�)�B
%
% [X,FVAL,EXITFLAG,OUTPUT] = FMINUNC(FUN,X0,...) �́A�J��Ԃ��� 
% OUTPUT.iterations�A�֐��̕]���� OUTPUT.funcCount�A�g�p�����A���S���Y�� 
% OUTPUT.algorithm�A(�g�p�����ꍇ) �������z�J��Ԃ��� OUTPUT.cgiterations�A
% (�g�p�����ꍇ) 1 ���̍œK�� OUTPUT.firstorderopt�A�I�����b�Z�[�W 
% OUTPUT.message �����\���� OUTPUT ���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINUNC(FUN,X0,...) �́A�� X �ł̊֐� 
% FUN �̌��z�l���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = FMINUNC(FUN,X0,...) �́A�� X 
% �ł̖ړI�֐� FUN �̃w�b�Z�s����o�͂��܂��B
% 
% ��
% FUN �́A@ ���g���Đݒ肷�邱�Ƃ��ł��܂��B
%
%        X = fminunc(@myfun,2)
%
% �����ŁAmyfun �́A���̂悤�ɕ\�킳��� MATLAB �֐��ł��B
% 
%        function F = myfun(x)
%         F = sin(x)+3;
%
% �^����ꂽ���z���g���āA�֐����ŏ������邽�߂ɁA���z�� 2 �Ԗڂ̈�����
% �Ȃ�悤�ɁA�֐� myfun ���C�����܂��B
% 
%        function [f,g] =  myfun(x)
%         f = sin(x) + 3;
%         g = cos(x);
% 
% �����āA(OPTIMSET ���g����) OPTIONS.GradObj �� 'on' �ɐݒ肵�A���z��
% �g�p�ł���悤�ɂ��܂��B
% 
%        options = optimset('GradObj','on');
%        x = fminunc(@myfun,4,options);
%
% FUN �́A�����֐��Ƃ��Ă��\���ł��܂��B
%
%        x = fminunc(@(x) 5*x(1)^2 + x(2)^2,[5;1])
%
% FUN ���p�����[�^�����ꂽ�ꍇ�A���Ɉˑ������p�����[�^���w�肵�Ė���
% �֐����g�p�ł��܂��B2 �Ԗڂ̈��� c �Ńp�����[�^�����ꂽ�֐� myfun ��
% �ŏ�������Ɖ��肵�܂��B�����ŁAmfun �͂��̂悤�� M-�t�@�C���֐��ł��B
%
%     function [f,g] = myfun(x,c)
%
%     f = c*x(1)^2 + 2*x(1)*x(2) + x(2)^2; % �֐�
%     g = [2*c*x(1) + 2*x(2)               % ���z
%          2*x(1) + 2*x(2)];
%
% c �̎w��l�ɑ΂��čœK�����邽�߂ɂ́A�ŏ��� c �ɐݒ肵�܂��B���ɁA
% x �������Ƃ��閳���֐� myfun ���`���A�ŏI�I�ɁAFMINUNC �ɓn���܂��B
%
%     c = 3;                              % �ŏ��Ƀp�����[�^���`���܂��B
%     options = optimset('GradObj','on'); % ���z���w�����܂��B
%     x = fminunc(@(x) myfun(x,c),[1;1],options)
%
%   �Q�l OPTIMSET, FMINSEARCH, FMINBND, FMINCON, @, INLINE.


%   Copyright 1990-2006 The MathWorks, Inc.
