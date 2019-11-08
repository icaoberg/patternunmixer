% FSOLVE   最小二乗法を使った多変数非線形方程式の求解
%
% FSOLVE は、つぎの形式の問題を解きます。
% 
% F(X) = 0    ここで、F と X は、ベクトルまたは行列です。  
%
% X = FSOLVE(FUN,X0) は、初期値を行列 X0 として、関数 FUN により表される
% 方程式を解きます。関数 FUN は、X を入力として、X で計算されるベクトル 
% (行列) の関数値 F を出力します。
%
% X = FSOLVE(FUN,X0,OPTIONS) は、オプション (OPTIONS) を設定して実行
% できます。ここで OPTIONS は OPTIMSET 関数で設定できる構造体です。詳細は、
% OPTIMSET を参照してください。Display, TolX, TolFun, DerivativeCheck, 
% Diagnostics, FunValCheck, Jacobian, JacobMult, JacobPattern, LineSearchType,
% NonlEqnAlgorithm, MaxFunEvals, MaxIter, OutputFcn, DiffMinChange, DiffMaxChange, 
% LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX のオプションが
% 使用できます。オプション Jacobian を使って、FUN を呼び出し、2 番目の
% 出力引数 J に、点 X でのヤコビ行列を設定できます。FUN は、X が長さ 
% n のときに、m 要素のベクトル F を出力する場合、J は、m 行 n 列の行列に
% なります。ここで、J(i,j) は、F(i) の x(j) による偏微分係数です 
% (ヤコビ行列 J は、F の勾配を転置したものです)。
%
% [X,FVAL] = FSOLVE(FUN,X0,...) は、X での方程式 FUN の値を出力します。
%
% [X,FVAL,EXITFLAG] = FSOLVE(FUN,X0,...) は、FSOLVE の終了状況を示す
% 文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する終了状況は、
% つぎの通りです。
%
%   1  FSOLVE は、解 X に収束したことを示します。
%   2  X の変化が指定した許容範囲より小さいことを示します。
%   3  目的関数値の変化が指定した許容範囲より小さいことを示します。
%   4  探索方向の大きさが指定した許容誤差より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを示します。
%  -1  最適化が、出力関数で終了していることを示します。
%  -2  アルゴリズムが解でない点に収束していることを示します。
%  -3  信頼領域半径が小さくなりすぎていることを示します。
%  -4  ライン探索は、現在の探索方向に沿って十分に残差を減少できないことを示します。
%
% [X,FVAL,EXITFLAG,OUTPUT] = FSOLVE(FUN,X0,...) は、繰り返し回数 
% OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、使用したアルゴリズム 
% OUTPUT.algorithm、(使用した場合) 共役勾配繰り返し回数 OUTPUT.cgiterations、
% (使用した場合) 1 次の最適性 OUTPUT.firstorderopt、終了メッセージ 
% OUTPUT.message をもつ構造体 OUTPUT を出力します。
%
% [X,FVAL,EXITFLAG,OUTPUT,JACOB] = FSOLVE(FUN,X0,...)  は、X での関数
% FUN のヤコビ行列を出力します。
% 
% 例
% FUN は、@ を使って設定することができます。
%
%        x = fsolve(@myfun,[2 3 4],optimset('Display','iter'))
%
% ここで、myfun は、つぎのように表される MATLAB 関数です。
%
%       function F = myfun(x)
%       F = sin(x);
%
% FUN は、無名関数としても表現できます。
%
%       x = fsolve(@(x) sin(3*x),[1 4],optimset('Display','off'))
%
% FUN がパラメータ化された場合、問題に依存したパラメータを指定して無名
% 関数を使用できます。2 番目の引数 c でパラメータ化された関数 myfun を
% 最小化すると仮定します。ここで、myfun はつぎのようなM-ファイル関数です。
%     
%       function F = myfun(x,c)
%       F = [ 2*x(1) - x(2) - exp(c*x(1))
%             -x(1) + 2*x(2) - exp(c*x(2))];
%
% c の指定値に対して最適化するためには、最初に c に設定します。次に、
% x を引数とする無名関数 myfun を定義し、最終的に、FSOLVE に渡します。
%
%       c = -1; % 最初にパラメータを定義します。
%       x = fsolve(@(x) myfun(x,c),[-5;-5])
%
%   参考 OPTIMSET, LSQNONLIN, @, INLINE.


%   Copyright 1990-2006 The MathWorks, Inc.
