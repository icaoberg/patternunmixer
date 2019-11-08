% QUADPROG  二次計画法
%
% X = QUADPROG(H,f,A,b) は、つぎの形式の二次計画問題を解きます。
%
% A*x <=  b のもとで、0.5*x'*H*x + f'*x を x に関して最小化します。
% 
% X = QUADPROG(H,f,A,b,Aeq,beq) は、線形等式制約 Aeq*x = beq のもとで、
% 上記の問題を解きます。
%
% X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB) は、設計変数 X の上下限の範囲を
% 与えることが可能です。この場合、LB <= X <= UB の範囲内で最適解が
% 探索されます。範囲の制約がない場合、LB と UB に空行列を設定してください。
% X(i) に下限がない場合、LB(i) = -Inf と設定し、X(i) に上限がない場合、
% UB(i) = Inf と設定してください。
%
% X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0) は、X0 を初期値とします。
%
% X = QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS)  は、オプション (OPTIONS) 
% を設定して実行できます。ここで OPTIONS は OPTIMSET 関数で設定できる構造体
% です。詳細は、OPTIMSET を参照してください。Display, Diagnostics, TolX, 
% TolFun, HessMult, LargeScale, MaxIter, PrecondBandWidth, TypicalX, 
% TolPCG, MaxPCGIter のオプションが使用できます。'final' と 'off' のみが
% Display に対して使われます ('iter' は使用できません)。
%
% X=QUADPROG(Hinfo,f,A,b,Aeq,beq,LB,UB,X0,OPTIONS,P1,P2,...) は、
% OPTIMSET('HessMult',HMFUN) が設定された場合、問題に依存するパラメータ 
% P1, P2,... を直接 HMFUN 関数に渡します。HMFUN は、ユーザによって与え
% られます。デフォルト値を使用する場合は、A, b, Aeq, beq, LB, UB, XO, 
% OPTIONS に空行列を渡してください。
%
% [X,FVAL] = QUADPROG(H,f,A,b) は、X での目的関数 FVAL = 0.5*X'*H*X + f'*X 
% の値を出力します。
%
% [X,FVAL,EXITFLAG] = QUADPROG(H,f,A,b) は、QUADPROG の終了状況を示す
% 文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する終了状況は、
% つぎの通りです。
%
%   1  QUADPROG は、解 X に収束したことを示します。
%   3  目的関数値の変化が指定した許容範囲より小さいことを示します。
%   4  局所最小値が見つかったことを示します。
%   0  繰り返し回数が最大回数に達していることを示します。
%  -2  可解が見つからなかったことを示します。
%  -3  問題に制約がないことを示します。
%  -4  条件が悪いため、最適化を中断したことを示します。
%  -7  探索方向の大きさが小さくなりすぎていることを示します。これ以上
%      計算を行なうことができません。
%
% [X,FVAL,EXITFLAG,OUTPUT] = QUADPROG(H,f,A,b) は、繰り返し回数 
% OUTPUT.iterations、使用したアルゴリズム OUTPUT.algorithm、(使用した
% 場合) 共役勾配繰り返し回数 OUTPUT.cgiterations 、(大規模法のみ) 1 次の
% 最適値 OUTPUT.firstorderopt、終了メッセージがある場合は OUTPUT.message を 
% もつ構造体 OUTPUT を出力します。
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = QUADPROG(H,f,A,b) は、最適解での
% ラグランジェ乗数 LAMBDA を出力します。 LAMBDA.ineqlin に
% 線形不等式 A を、LAMBDA.eqlin に線形等式 Aeq を、LAMBDA.lower に 
% LB を、LAMBDA.upper に UB を設定しています。


%   Copyright 1990-2006 The MathWorks, Inc. 
