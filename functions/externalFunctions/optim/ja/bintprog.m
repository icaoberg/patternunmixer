% BINTPROG  2 �l�����v����
%
%   BINTPROG �́A2 �l�����v����������܂��B
%   ���̐���̂��ƂŁAf'*X �� X �Ɋւ��čŏ������܂��B
%                    A*X <= b,
%                    Aeq*X = beq,
%                    �����ŁAX �̗v�f�� 2 �l�����A
%                    ���Ȃ킿 0 �� 1 �ł��B
%
%   X = BINTPROG(f) �́A f'*X ���ŏ���������������܂��B�����ŁAX ��
%   �v�f�́A2 �l�����ł��B
%
%   X = BINTPROG(f,A,b) �́A���`�s���� A*X <= b �̂��ƂŁAf'*X ���ŏ���
%   ������������܂��B�����ŁAX �̗v�f�́A2 �l�����ł��B
%
%   X = BINTPROG(f,A,b,Aeq,beq) �́A���`���� Aeq*X = beq �� ���`�s���� 
%   A*X <= b �̂��Ƃ� f'*X ���ŏ���������������܂��B�����ŁAX �̗v�f�́A
%   2 �l�����ł��B
%
%   X = BINTPROG(f,A,b,Aeq,beq,X0) �́A�����l�� X0 �ɐݒ肵�܂��B�����l 
%   X0 �́A2 �l�����ŁA�����łȂ���΂Ȃ�܂���B�����łȂ��ꍇ�A
%   ��������܂��B
%
%   X = BINTPROG(f,A,b,Aeq,beq,X0,OPTIONS) �́A�I�v�V���� (OPTIONS) ��
%   �ݒ肵�Ď��s�ł��܂��B�����ŁAOPTIONS �� OPTIMSET �֐��Őݒ�ł���
%   �\���̂ł��B�ڍׂ́AOPTIMSET ���Q�Ƃ��Ă��������BBranchStrategy, 
%   Diagnostics, Display, NodeDisplayInterval, MaxIter, MaxNodes, 
%   MaxRLPIter, MaxTime, NodeSearchStrategy, TolFun, TolXInteger, 
%   TolRLPFun �̃I�v�V�������g�p�ł��܂��B
%
%   [X,FVAL] = BINTPROG(...) �́A�� X �ł̖ړI�֐��̒l FVAL = f'*X ��
%   �o�͂��܂��B
%
%   [X,FVAL,EXITFLAG] = BINTPROG(...) �́ABINTPROG �̏I���󋵂����������� 
%   EXITFLAG ��Ԃ��܂��BEXITFLAG �̉\�Ȓl�ƑΉ�����I���󋵂́A����
%   �ʂ�ł��B
%
%      1  BINTPROG �́A�� X �Ɏ����������Ƃ������܂��B
%      0  �J��Ԃ��񐔂��ő�񐔂ɒB���Ă��邱�Ƃ������܂��B
%     -2  ����������Ȃ��������Ƃ������܂��B
%     -4  �������邱�ƂȂ� MaxNodes �ɒB�������Ƃ������܂��B
%     -5  �������邱�ƂȂ� MaxTime �ɒB�������Ƃ������܂��B
%     -6  LP �ɘa�����������߂̃m�[�h�Ŏ��s���ꂽ�J��Ԃ��񐔂��A
%         �������邱�ƂȂ� MaxRLPIter �ɒB�������Ƃ������܂��B
%
%   [X,FVAL,EXITFLAG,OUTPUT] = BINTPROG(...) �́A�J��Ԃ��� OUTPUT.iterations�A
%   �T�������m�[�h�̐� OUTPUT.nodes�A���s���� (�b�P��) OUTPUT.time�A
%   �g�p�����A���S���Y�� OUTPUT.algorithm�A����@ OUTPUT.branchStrategy�A
%   �m�[�h�������@ OUTPUT.nodeSrchStrategy�A�I�����b�Z�[�W OUTPUT.message 
%   �����\���� OUTPUT ���o�͂��܂��B
%
%   ��
%     f = [-9; -5; -6; -4];
%     A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
%     b = [9; 1; 0; 0];
%     X = bintprog(f,A,b)


%   Copyright 1990-2006 The MathWorks, Inc.
