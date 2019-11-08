% FGOALATTAIN   多目的ゴール到達最適化問題
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT) は、X を変化させることにより、関数 
% FUN で設定される目的関数 (F) を、ゴール（GOAL）に到達させます。ゴールは、
% WEIGHT に従って、重み付けられます。これを行なうときに、つぎの形式の
% 非線形計画問題を解きます。
% 
% LAMBDA :  F(X)-WEIGHT.*LAMBDA<=GOAL を X, LAMBDA に関して最小化します。
%           
% 関数 FUN は、X を入力として、X で計算されるベクトル (行列) の目的値 
% F を出力します。X0 は、スカラ、ベクトル、または行列です。
% 
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B) は、線形不等式制約 A*X <= B の
% もとで、ゴール到達問題を解きます。
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq) は、線形等式制約 
% Aeq*X = Beq も考慮して、ゴール到達問題を解きます。
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB)  は、設計変数 X 
% の上下限を与えることが可能です。この場合は、LB <= X <= UB の範囲内で
% 最適解が探索されます。範囲の設定を行なわない場合は、LB と UB に空行列を
% 設定してください。また、 X(i) に下限がない場合、LB(i) = -Inf とし、
% X(i) に上限がない場合、UB(i) = Inf と設定してください。
% 
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON) は、
% NONLCON (通常は、M-ファイル NONLCON.m) で定義した制約のもとで、ゴール
% 到達問題を解きます。ここで関数 NONLCON は、
% feval: [C, Ceq] = feval(NONLCON,X) のように呼び出される場合、それぞれ、
% 非線形の不等式制約と等式制約を表わす C と Ceq ベクトルを返します。
% FGOALATTAIN は、C(X)< = 0 と Ceq(X) = 0 になるように最適化します。
%
% X = FGOALATTAIN(FUN,X0,GOAL,WEIGHT,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) 
% は、オプション (OPTIONS) を設定して実行できます。ここで、OPTIONS は、
% OPTIMSET 関数で設定できる構造体です。詳細は、OPTIMSET を参照してください。
% Display, TolX, TolFun, TolCon, DerivativeCheck, FunValCheck, GradObj, 
% GradConstr, MaxFunEvals, MaxIter, MeritFunction, GoalsExactAchieve, 
% Diagnostics, DiffMinChange, DiffMaxChange, OutputFcn, TypicalX の
% オプションを使用できます。GradObj オプションを使って、FUN を呼び出し、
% 2 番目の出力引数 G に X での偏微分係数 df/dX を設定できます。
% [F,G] = feval(FUN,X) の形式で呼び出します。GradConstr オプションを使って、
% つぎのように、4 つの出力引数で呼び出される NONLCON を指定できます。
% [C,Ceq,GC,GCeq] = feval(NONLCON,X) の形式で呼び出します。ここで、GC は、
% 不等号制約ベクトル C の偏微分係数、GCeq は、等式制約ベクトル Ceq の
% 偏微分係数です。オプションを設定しない場合は、OPTIONS = [] を使用して
% ください。
%
% [X,FVAL] = FGOALATTAIN(FUN,X0,...) は、解 X での目的関数 FUN の値を
% 出力します。
%
% [X,FVAL,ATTAINFACTOR] = FGOALATTAIN(FUN,X0,...) は、解 X での到達ファクタ
% を出力します。ATTAINFACTOR が負の場合、ゴールは過到達になります。また、
% 正の場合、欠到達になります。
%
% [X,FVAL,ATTAINFACTOR,EXITFLAG] = FGOALATTAIN(FUN,X0,...) は、FGOALATTAIN の
% 終了状況を示す文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する
% 終了状況は、つぎの通りです。
%
%   1  FGOALATTAIN は、解 X に収束していることを示します。
%   4  探索方向の大きさが指定した許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   5  探索方向に沿った関数の勾配の大きさが許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを示します。
%  -1  最適化が、出力関数で終了していることを示します。
%  -2  可解が見つからなかったことを示します。
%   
% [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT] = FGOALATTAIN(FUN,X0,...) は、
% 繰り返し数 OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、
% 最終ステップのノルム OUTPUT.stepsize、最終のライン探索のステップ長 
% OUTPUT.lssteplength、使用したアルゴリズム OUTPUT.algorithm、1 次の最適性 
% OUTPUT.firstorderopt、終了メッセージ  OUTPUT.message をもつ構造体 
% OUTPUT を出力します。
% 
% [X,FVAL,ATTAINFACTOR,EXITFLAG,OUTPUT,LAMBDA] = FGOALATTAIN(FUN,X0,...)
% は、解でのラグランジェ乗数 LAMBDA を出力します。LAMBDA 構造体は、
% LAMBDA.lower に LB を、LAMBDA.upper に UB を、LAMBDA.ineqlin に線形
% 不等式を、LAMBDA.eqlin に線形等式を、LAMBDA.ineqnonlin に非線形不等式を、
% LAMBDA.eqnonlin に非線形等式を設定します。
%
% 詳細は、M-ファイル FGOALATTAIN.M を参照してください。
%
% 参考 OPTIMSET, OPTIMGET.


%   Copyright 1990-2006 The MathWorks, Inc.
