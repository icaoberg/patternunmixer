% LSQLIN   制約付き線形最小二乗問題
%
% X = LSQLIN(C,d,A,b) は、つぎの形式の不等号の制約条件付き最小二乗問題を
% 解きます。
%
% A*x <=  b の制約のもとで、0.5*(NORM(C*x-d)).^2 を x に関して最小化します。
% ここで、C は、m 行 n 列の行列です。
%
% X = LSQLIN(C,d,A,b,Aeq,beq) は、つぎの形式の (等号の制約付き) 最小二乗
% 問題を解きます。
%
% A*x <=  b と Aeq*X = Beq を考慮して、0.5*(NORM(C*x-d)).^2 を x に関して、
% 最小化します。
%
% X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB) は、設計変数 X の上下限の範囲を与える
% ことが可能です。この場合、LB <= X <= UB の範囲内で最適解が探索されます。
% 範囲の制約がない場合、LB と UB に空行列を設定してください。X(i) に
% 下限がない場合、LB(i) = -Inf と設定し、X(i) に上限がない場合、UB(i) = Inf 
% と設定してください。
%
% X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0) は、X0 を初期値とします。
%
% X = LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0,OPTIONS) は、オプション (OPTIONS) 
% を設定して実行できます。ここで OPTIONS は OPTIMSET 関数で設定できる
% 構造体です。詳細は、OPTIMSET を参照してください。Display, Diagnostics, 
% TolFun, LargeScale, MaxIter, JacobMult, PrecondBandWidth, TypicalX, 
% TolPCG, MaxPCGIter のオプションを使用できます。'final' と 'off' のみが、
% Display に対して使われます ('iter' は使用できません)。
%
% X=LSQLIN(C,d,A,b,Aeq,beq,LB,UB,X0,OPTIONS,P1,P2,...) は、
% OPTIMSET('JacobMult',JMFUN) が設定された場合、問題に依存するパラメータ 
% P1,P2,... を直接 JMFUN 関数に渡します。JMFUN は、ユーザによって指定
% します。デフォルト値を使用するには、A, b, Aeq, beq, LB, UB, XO, 
% OPTIONS に空行列を渡してください。
%
% [X,RESNORM] = LSQLIN(C,d,A,b) は、残差の二乗 2 ノルム値 norm(C*X-d)^2 を
% 出力します。
%
% [X,RESNORM,RESIDUAL] = LSQLIN(C,d,A,b) は、残差 C*X-d を出力します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG] = LSQLIN(C,d,A,b) は、LSQLIN の終了状況を
% 示す文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する終了状況は、
% つぎの通りです。
%
%   1  LSQLIN は、解 X に収束したことを示します。
%   3  残差の変化が指定した許容範囲より小さいことを示します。
%   0  繰り返し回数が最大回数に達していることを示します。
%  -2  問題が不可解であることを示します。
%  -4  条件が悪いため、最適化を中断したことを示します。
%  -7  探索方向の大きさが小さくなりすぎていることを示します。これ以上
%      計算を行なうことができません。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQLIN(C,d,A,b) は、繰り返し 
% OUTPUT.iterations、使用したアルゴリズム OUTPUT.algorithm、(使用した場合) 
% 共役勾配繰り返し回数 OUTPUT.cgiterations、(大規模問題のみ) 1 次の最適性 
% OUTPUT.firstorderopt、終了メッセージ OUTPUT.message をもつ構造体 OUTPUT を
% 出力します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQLIN(C,d,A,b)  は、
% 最適解でのラグランジェ乗数 LAMBDA を出力します。LAMBDA.ineqyalities に線形
% 不等式 C を、LAMBDA.eqlin に線形等式 Ceq を、LAMBDA.lower に LB を、
% LAMBDA.upper に UB を設定します。


%   Copyright 1990-2006 The MathWorks, Inc.
