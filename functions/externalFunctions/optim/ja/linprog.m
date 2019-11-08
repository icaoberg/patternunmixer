% LINPROG  ���`�v����
%
% X = LINPROG(f,A,b) �́A���̌`���̐��`�v����������܂��B
% A*x <= b �̂��ƂŁAf'*x �� x �Ɋւ��čŏ������܂��B
%  
% X = LINPROG(f,A,b,Aeq,beq) �́A�������� Aeq*X = Beq ���l�����āA��L��
% ���������܂��B
% 
% X = LINPROG(f,A,b,Aeq,beq,LB,UB) �́A�݌v�ϐ� X �̏㉺���͈̔͂�^����
% ���Ƃ��\�ł��B���̏ꍇ�ALB <= X <= UB �͈͓̔��ōœK�����T������܂��B
% �͈͂̐��񂪂Ȃ��ꍇ�ALB �� UB �ɋ�s���ݒ肵�Ă��������BX(i) ��
% �������Ȃ��ꍇ�ALB(i) = -Inf �Ɛݒ肵�AX(i) �ɏ�����Ȃ��ꍇ�AUB(i) = Inf 
% �Ɛݒ肵�Ă��������B
%
% X = LINPROG(f,A,b,Aeq,beq,LB,UB,X0) �́A�����l�� X0 �Ƃ��܂��B����
% �I�v�V�����́Aactive-set �A���S���Y�����g�p����Ƃ��̂ݗ��p�ł��܂��B
% �f�t�H���g�̓��_�A���S���Y���́A�S�Ă̋�łȂ������l�𖳎����܂��B
%
% X=LINPROG(f,A,b,Aeq,beq,LB,UB,X0,OPTIONS) �́A�I�v�V���� (OPTIONS) 
% ��ݒ肵�Ď��s�ł��܂��B������ OPTIONS �� OPTIMSET �֐��Őݒ�ł���
% �\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BDisplay, Diagnostics, 
% TolFun, LargeScale, MaxIter �̃I�v�V�������g�p�ł��܂��BLargeScale �� 
% 'off' �̏ꍇ�A'final' �� 'off' �݂̂��ADisplay �ɑ΂��Ďg���܂� 
% (LargeScale �� 'on' �̏ꍇ�A'iter' ���g�p�ł��܂�)�B
%
% [X,FVAL] = LINPROG(f,A,b) �́A�� X �ł̖ړI�֐��̒l FVAL = f'*X ��
% �o�͂��܂��B
%
% [X,FVAL,EXITFLAG] = LINPROG(f,A,b) �́ALINPROG �̏I���󋵂����������� 
% EXITFLAG ���o�͂��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A����
% �ʂ�ł��B
%
%   1  LINPROG �́A�� X �Ɏ����������Ƃ������܂��B
%   0  �J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%  -2  ����������Ȃ��������Ƃ������܂��B
%  -3  ���ɐ��񂪂Ȃ����Ƃ������܂��B
%  -4  �A���S���Y���̕]������ NaN �̒l�����������Ƃ������܂��B
%  -5  �o�΁A����̗������s���ł��邱�Ƃ������܂��B
%  -7  �T�������̑傫�����������Ȃ肷���Ă��邱�Ƃ������܂��B����ȏ�
%      �v�Z���s�Ȃ����Ƃ��ł��܂���B
%
% [X,FVAL,EXITFLAG,OUTPUT] = LINPROG(f,A,b) �́A�J��Ԃ��� OUTPUT.iterations�A
% �g�p�����A���S���Y�� OUTPUT.algorithm�A(�g�p�����ꍇ) �������z�J��Ԃ��� 
% OUTPUT.cgiterations�A�I�����b�Z�[�W OUTPUT.message �����\���� OUTPUT 
% ���o�͂��܂��B
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = LINPROG(f,A,b) �́A�� X �ł̃��O�����W�F
% �搔 LAMBDA ���o�͂��܂��BLAMBDA �\���̂́ALAMBDA.lower �� LB ���A
% LAMBDA.upper �� UB ���ALAMBDA.ineqlin �ɐ��`�s�������ALAMBDA.eqlin ��
% ���`�������ALAMBDA.ineqnonlin �ɔ���`�s�������ALAMBDA.eqnonlin ��
% ����`������ݒ肵�܂��B
%   
% ���� �F LINPROG �̑�K�͖�� (�f�t�H���g) �́A��o�Ζ@���g�p���܂��B
% ����Ƒo�Ζ��́A���Ɏ����ɑ΂��ĕs���łȂ���΂Ȃ�܂���B����A
% �o�Ζ��A���邢�͗����̂����ꂩ�����ł���Ƃ������b�Z�[�W�́A�K�؂�
% �^�����܂��B����̕W���I�Ȍ`���́A���̒ʂ�ł��B
% 
%      A*x = b, x >= 0 �̂��ƂŁAf'*x ���ŏ������܂��B
% 
% �o�Ζ��́A���̌`���ɂȂ�܂��B
% 
%      A'*y + s = f, s >= 0 �̂��ƂŁAb'*y ���ő剻���܂��B


%   Copyright 1990-2006 The MathWorks, Inc.
