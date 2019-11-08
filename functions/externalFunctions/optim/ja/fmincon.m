% FMINCON  制約付き多変数関数の最小化
%
% FMINCON は、つぎの形式の問題を解きます。
% つぎの制約のもとで、F(X) を、X に関して最小化します。
% 　　　　　　A*X  <= B, Aeq*X  = Beq (線形制約)
%       　　　C(X) <= 0, Ceq(X) = 0   (非線形制約)
%       　　　LB <= X <= UB            
% 
% X = FMINCON(FUN,X0,A,B) は、初期 X0 値で、線形不等式制約 A*X <= B の
% もとで、関数 FUN を最小化する X を求めます。FUN は、X を入力として、
% X で計算されるスカラの関数値を出力します。X0 は、スカラ、ベクトル、
% または、行列です。
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq) は、A*X <= B と同様に線形等式制約 
% Aeq*X = Beq を考慮して FUN を最小化します (不等式制約がない場合は、
% A=[], B=[] と設定します)。
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB) は、設計変数 X の上下限の範囲を
% 与えることが可能です。この場合、LB <= X <= UB の範囲内で最適解が探索
% されます。範囲の制約がない場合、LB と UB に空行列を設定してください。
% X(i) に下限がない場合、LB(i) = -Inf と設定し、X(i) に上限がない場合、
% UB(i) = Inf と設定してください。
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) は、NONLCON で定義された
% 制約のもとで最小化します。ここで関数 NONLCON は、X を入力として、
% 非線形不等式と非線形等式を表わす C と Ceq ベクトルを返します。FMINCON は、
% C(X) <= 0 と Ceq(X) = 0 となる  FUN を最小化します。(範囲の制約がない場合、
% LB=[] と UB=[] として設定してください。
%
% X = FMINCON(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS) は、オプション 
% (OPTIONS) を設定して実行できます。ここで、OPTIONS は OPTIMSET 関数で
% 設定できる構造体です。詳細は、OPTIMSET を参照してください。
% Display, TolX, TolFun, TolCon, DerivativeCheck, Diagnostics, FunValCheck, 
% GradObj, GradConstr, Hessian, MaxFunEvals, MaxIter, DiffMinChange, 
% DiffMaxChange, LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX, 
% Hessian, HessMult, HessPattern, OutputFcn オプションが使用できます。GradObj 
% オプションを使って、2番目の出力引数 G に点 X での関数の偏微分係数 df/dX を
% 出力する FUN を設定します。オプション Hessian を使って、FUN を呼び出し、
% 3番目の引数 H に、点 X でのヘッセ行列 (2階微分) を設定できます。Hessian は、
% 大規模問題のみに使用され、ライン探索法では使用されません。GradConstr 
% オプションを使って、NONLCON が、3 番目の出力引数に GC 、4 番目の引数に 
% GCeq をもつように設定します。ここで、GC は不等号制約ベクトル C の偏微分
% 係数、GCeq は、等式制約ベクトル Ceq の偏微分係数です。オプションを設定
% しない場合は、OPTIONS = [] を使用してください。
% 
% [X,FVAL] = FMINCON(FUN,X0,...) は、解 X での目的関数 FUN の値を出力します。
%
% [X,FVAL,EXITFLAG] = FMINCON(FUN,X0,...) は、FMINCON の終了状況を示す
% 文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する終了状況は、
% つぎの通りです。
%
%   1  1 次の最適条件が指定した許容範囲を満足していることを示します。
%   2  X の変化が指定した許容範囲より小さいことを示します。
%   3  目的関数値の変化が指定した許容範囲より小さいことを示します。
%   4  探索方向の大きさが指定した許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   5  探索方向に沿った関数の勾配の大きさが許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを
%      示します。
%  -1  最適化が、出力関数で終了していることを示します。
%  -2  可解が見つからなかったことを示します。
% 
% [X,FVAL,EXITFLAG,OUTPUT] = FMINCON(FUN,X0,...) は、
% 繰り返し数 OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、最終
% ステップのノルム OUTPUT.stepsize、使用したアルゴリズム OUTPUT.algorithm、
% 1 次の最適性 OUTPUT.firstorderopt、終了メッセージ OUTPUT.message をもつ
% 構造体 OUTPUT を出力します。中規模アルゴリズムは、最終のライン探索の
% ステップ長 OUTPUT.lssteplength を返し、大規模アルゴリズムは、共役勾配
% 繰り返し回数 OUTPUT.cgiterations を出力します。
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FMINCON(FUN,X0,...)  は、解での 
% ラグランジェ乗数 LAMBDA を出力します。LAMBDA 構造体は、LAMBDA.lower に 
% LB を、LAMBDA.upper に UB を、LAMBDA.ineqlin に線形不等式を、LAMBDA.eqlin 
% に線形等式を、LAMBDA.ineqnonlin に非線形不等式を、LAMBDA.eqnonlin に
% 非線形等式を設定します。
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD] = FMINCON(FUN,X0,...) は、解 X で
% の関数 FUN の勾配値も出力します。
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = FMINCON(FUN,X0,...) は
% 解 X での関数 FUN のヘッセ行列を出力します。
%　
% 例：
% FUN は、@ を使って設定することができます : X = fmincon(@hump,...) 
% この場合、F = humps(X) は、X での HUMPS 関数のスカラ関数値 F を出力します。
% 
% FUN は、無名関数としても表現できます。
% 
%        X = fmincon(@(x) 3*sin(x(1))+exp(x(2)),[1;1],[],[],[],[],[0 0])
% 
% は、X = [0;0] を出力します。
% 
% FUN または NONLCON がパラメータ化された場合、問題に依存したパラメータを
% 指定して無名関数を使用できます。関数 myfun を最小化すると仮定し、また
% 非線形制約条件 mycon を考慮します。ここで、これらの 2 つの関数は 
% 2 番目の引数 a1 と a2 でそれぞれパラメータ化されます。ここで、myfun と 
% mycon はつぎのような M-ファイル関数です。
%
%        function f = myfun(x,a1)
%        f = x(1)^2 + a1*x(2)^2;
%
% そして、
%
%        function [c,ceq] = mycon(x,a2)
%        c = a2/x(1) - x(2);
%        ceq = [];
%
% 指定された a1 と a2 に対して最適化を行なうには、まず、これらのパラメータを
% 設定します。そして、a1 と a2 を 1 つの引数とする無名関数を 2 つ定義します。
% 最終的に、これらの無名関数を FMINCON に渡します。
%
%        a1 = 2; a2 = 1.5; % 最初にパラメータを定義します。
%        options = optimset('LargeScale','off'); % 中規模アルゴリズムを実行させます。
%        x = fmincon(@(x)myfun(x,a1),[1;2],[],[],[],[],[],[],@(x)mycon(x,a2),options)
%
%   参考 OPTIMSET, FMINUNC, FMINBND, FMINSEARCH, @, FUNCTION_HANDLE.


%   Copyright 1990-2006 The MathWorks, Inc.
