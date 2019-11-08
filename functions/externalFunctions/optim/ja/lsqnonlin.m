% LSQNONLIN   非線形最小二乗問題
%
% LSQNONLIN は、つぎの形式の問題を解きます。
% sum {FUN(X).^2} を、X に関して最小化します。ここで、X と関数 FUN で
% 返される値は、ベクトル、または、行列です。
%
% X = LSQNONLIN(FUN,X0) は、X0 を初期値として、関数 FUN に記述される関数
% の二乗和を最小化する X を求めます。関数 FUN は、X を入力として、X で
% 計算されるベクトル (または行列) の関数値 F を出力します。注意 ： FUN は、
% FUN(X) を返し、二乗和 sum(FUN(X).^2) は返しません (この処理は、アルゴリズム
% 内で陰的に行なわれます)。
%
% X = LSQNONLIN(FUN,X0,LB,UB) は、設計変数 X の上下限の範囲を与えることが
% 可能です。この場合、LB <= X <= UB の範囲で最適解が探索されます。範囲の
% 制約がない場合、LB と UB に空行列を設定してください。X(i) に下限がない場合、
% LB(i) = -Inf と設定し、X(i) に上限がない場合、UB(i) = Inf と設定してください。
% 
% X = LSQNONLIN(FUN,X0,LB,UB,OPTIONS) は、オプション (OPTIONS) を設定して
% 実行できます。ここで OPTIONS は OPTIMSET 関数で設定できる構造体です。
% 詳細は、OPTIMSET を参照してください。Display, TolX, TolFun, DerivativeCheck, 
% Diagnostics, FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType, 
% LevenbergMarquardt, MaxFunEvals, MaxIter, DiffMinChange, DiffMaxChange, 
% LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX, OutputFcn 
% のオプションが使用できます。オプション Jacobian を使って、FUN を呼び出し、
% 2 番目の出力引数 J に、点 X でのヤコビ行列を設定できます。X が長さ n の
% ときに、FUN が m 要素のベクトル F を返す場合、J は、m 行 n 列の行列に
% なります。ここで、J(i,j) は、F(i) の x(j) による偏微分係数です (ヤコビ行列 
% J は、F の勾配を転置したものです)。
%
% [X,RESNORM] = LSQNONLIN(FUN,X0,...) は、X での残差の 2 ノルム sum(FUN(X).^2) 
% を出力します。
%
% [X,RESNORM,RESIDUAL] = LSQNONLIN(FUN,X0,...) は、解 X での残差値 
% RESIDUAL = FUN(X) を返します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG] = LSQNONLIN(FUN,X0,...) は、LSQNONLIN の
% 終了状況を示す文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する
% 終了状況は、つぎの通りです。
%
%   1  LSQNONLIN は、解 X に収束したことを示します。
%   2  X の変化が指定した許容範囲より小さいことを示します。
%   3  目的関数値の変化が指定した許容範囲より小さいことを示します。
%   4  探索方向の大きさが指定した許容誤差より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを示します。
%  -1  アルゴリズムが、出力関数で終了していることを示します。
%  -2  境界が矛盾していることを示します。
%  -4  ライン探索は、現在の探索方向に沿って十分に残差を減少できないことを示します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQNONLIN(FUN,X0,...) は、繰り
% 返し回数 OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、使用した
% アルゴリズム OUTPUT.algorithm、(使用した場合) 共役勾配繰り返し回数 
% OUTPUT.cgiterations、(使用した場合) 1 次の最適性 OUTPUT.firstorderopt、
% 終了メッセージ OUTPUT.message をもつ構造体 OUTPUT を出力します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQNONLIN(FUN,X0,...) は、
% 解でのラグランジェ乗数 LAMBDA を出力します。LAMBDA.lower に LB を、
% LAMBDA.upper にUB を設定します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = LSQNONLIN(FUN,X0,...)
% は、X での関数 FUN のヤコビ行列を出力します。
%
% 例
% FUN は、@ を使って設定することができます。
%        x = lsqnonlin(@myfun,[2 3 4])
%
% ここで、myfun は、つぎのように表される MATLAB 関数です。
%
%       function F = myfun(x)
%       F = sin(x);
%
% FUN は、無名関数としても表現できます。
%
%       x = lsqnonlin(@(x) sin(3*x),[1 4])
%
% FUN または NONLCON がパラメータ化された場合問題に依存したパラメータを
% 指定して無名関数を使用できます。2 番目の引数 c でパラメータ化された関数 
% myfun を最小化すると仮定します。ここで、mfun はつぎのような M-ファイル
% 関数です。
%
%       function F = myfun(x,c)
%       F = [ 2*x(1) - exp(c*x(1))
%             -x(1) - exp(c*x(2))
%             x(1) - x(2) ];
%
% c の指定値に対して最適化するためには、最初に c に設定します。次に、
% x を引数とする無名関数 myfun を定義し、最終的に、LSQNONLIN に渡します。
%
%       c = -1; % 最初にパラメータを定義します。
%       x = lsqnonlin(@(x) myfun(x,c),[1;1])
%
%   参考 OPTIMSET, LSQCURVEFIT, FSOLVE, @, INLINE.


%   Copyright 1990-2006 The MathWorks, Inc.
