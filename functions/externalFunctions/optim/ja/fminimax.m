% FMINIMAX   多変数関数のミニマックス問題
%
% FMINIMAX は、つぎの形式の問題を解きます。
% {FUN(X)} を X に関して最小化します。ここで、FUN と X は、ベクトル、または
% 行列です。
%
% X = FMINIMAX(FUN,X0) は、X0 を初期値として、関数 FUN のミニマックス解 
% X を求めます。FUN は、X を入力として、X で計算されるベクトル (行列) の
% 関数値 F を出力します。X0 は、スカラ、ベクトル、または、行列です。
%
% X = FMINIMAX(FUN,X0,A,B) は、線形不等式制約 A*X <= B のもとでミニマックス
% 問題を解きます。
%
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq) は、線形等式制約 Aeq*X = Beq を考慮して
% ミニマックス問題を解きます (不等式制約がない場合には、A=[],B=[] と設定
% します)。
% 
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB) は、設計変数 X の上下限の範囲を
% 与えることが可能です。ここで、LB <= X <= UB 範囲内で最適解が探索されます。
% 範囲の制約がない場合、LB と UB に空行列を設定してください。X(i) に下限が
% ない場合、LB(i) = -Inf と設定し、X(i) に上限がない場合、UB(i) = Inf と
% 設定してください。
% 
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) は、NONLCON (通常は 
% M-ファイル NONLINCON.m) で定義された制約のもとで、ミニマックス問題を
% 解きます。ここで関数 NONLCON は、feval: [C, Ceq] = feval(NONLCON,X) の
% ように呼び出される場合、それぞれ、非線形の不等式制約と等式制約を表わす 
% C と Ceq ベクトルを返します。FGOALATTAIN は、C(X)< = 0 と Ceq(X) = 0 に
% なるように最適化します。
%
% X = FMINIMAX(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)  は、オプション 
% (OPTIONS) を設定して実行できます。ここで、OPTIONS は OPTIMSET 関数で
% 設定できる構造体です。詳細は、OPTIMSET を参照してください。Display, 
% TolX, TolFun, TolCon, DerivativeCheck, FunValCheck, GradObj, GradConstr, 
% MaxFunEvals, MaxIter, MeritFunction, MinAbsMax, Diagnostics, DiffMinChange, 
% DiffMaxChange, OutputFcn, TypicalX のオプションが使用できます。GradObj 
% オプションを使って、FUN を呼び出し、2 番目の出力引数 G に点 X での
% 偏微分係数 df/dX を設定できます。[F,G] = feval(FUN,X) の形式で呼び出します。
% GradConstr オプションを使って、つぎのように、4つの出力引数で呼び出される 
% NONLCON を指定できます。[C,Ceq,G,C,GCeq] = feval(NONLINCON,X) の形式で
% 呼び出します。ここで、GC は、不等号制約ベクトル C の偏微分係数、GCeq は、
% 等式制約ベクトル Ceq の偏微分係数です。オプションを設定しない場合は、
% OPTIONS = [] を使用してください。
%
% [X,FVAL] = FMINIMAX(FUN,X0,...) は、解 X での目的関数 FVAL = feval(FUN,X) 
% の値を出力します。
%
% [X,FVAL,MAXFVAL] = FMINIMAX(FUN,X0,...) は、解 X での MAXFVAL = max { FUN(X) } 
% を出力します。
%
% [X,FVAL,MAXFVAL,EXITFLAG] = FMINIMAX(FUN,X0,...) は、FMINIMAX の終了状況を
% 示す文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する終了状況は、
% つぎの通りです。
%
%   1  1 次の最適条件が指定した許容範囲を満足していることを示します。
%   4  探索方向の大きさが指定した許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   5  探索方向に沿った関数の勾配の大きさが許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを示します。
%  -1  最適化が、出力関数で終了していることを示します。
%  -2  可解が見つからなかったことを示します。
%   
% [X,FVAL,MAXFVAL,EXITFLAG,OUTPUT] = FMINIMAX(FUN,X0,...) は、実行した
% 繰り返し数 OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、最終
% ステップのノルム OUTPUT.stepsize、最終のライン探索のステップ長 
% OUTPUT.lssteplength、使用したアルゴリズム OUTPUT.algorithm、1 次の最適性 
% OUTPUT.firstorderopt、終了メッセージ  OUTPUT.message をもつ構造体 
% OUTPUT を出力します。
%
% [X,FVAL,MAXFVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINIMAX(FUN,X0,...)  は、
% 解でのラグランジェ乗数 LAMBDA を出力します。LAMBDA 構造体は、
% LAMBDA.lower に LB を、LAMBDA.upper に UB を、LAMBDA.ineqlin に線形
% 不等式を、LAMBDA.eqlin に線形等式を、LAMBDA.ineqnonlin に非線形不等式を、
% LAMBDA.eqnonlin に非線形等式を設定します。
%
% 例
% FUN は、@ を使って設定することができます。
%
%        x = fminimax(@myfun,[2 3 4])
%
% ここで、myfun は、つぎのように表される MATLAB 関数です。
% 
%   function F = myfun(x)
%   F = cos(x);
% 
% FUN は、無名関数としても表現できます。
% 
%       x = fminimax(@(x) sin(3*x),[2 5])
%
% FUN がパラメータ化された場合、問題に依存したパラメータを指定して無名
% 関数を使用できます。2 番目の引数 c でパラメータ化された関数 myfun を
% 最小化すると仮定します。ここで、mfun はつぎのような M-ファイル関数です。
%
%       function F = myfun(x,c)
%       F = [x(1)^2 + c*x(2)^2;
%            x(2) - x(1)];
%
% c の指定値に対して最適化するためには、最初に c に設定します。次に、
% x を引数とする無名関数 myfun を定義し、最終的に、FMINIMAX に渡します。
%
%       c = 2; % 最初にパラメータを定義します。
%       x = fminimax(@(x) myfun(x,c),[1;1])
%
%   参考 OPTIMSET, @, INLINE, FGOALATTAIN, LSQNONLIN.


%   Copyright 1990-2006 The MathWorks, Inc.
