% LSQCURVEFIT   非線形最小二乗問題 (回帰分析に特化)
%
% LSQCURVEFIT は、つぎの形式の問題を解きます。
% sum {(FUN(X,XDATA)-YDATA).^2} を X に関して最小化します。
%                               ここで、X, XDATA, YDATA, FUN の
%                               出力値は、ベクトルまたは行列です。
%
% X = LSQCURVEFIT(FUN,X0,XDATA,YDATA) は、X0 を初期値とし、関数 FUN で
% 定義する非線形関数をデータ YDATA に、最小二乗的に最適近似する係数 X を求め
% ます。FUN は、X と XDATA を入力として、ベクトル (または、行列) 関数値 
% F を出力します。ここで、F は、YDATA と同じサイズで、X と XDATA で計算した
% ものです。注意：FUN は、FUN(X,XDATA) を返し、sum(FUN(X,XDATA).^2) で
% 表される二乗和ではありません。この処理は、アルゴリズムの中で、陰的に
% 行なわれます。
%
% X = LSQCURVEFIT(FUN,X0,XDATA,YDATA,LB,UB) は、設計変数 X の上下限の
% 範囲を与えることが可能です。ここで、LB <= X <= UB の範囲内で最適解が
% 探索されます。範囲の制約がない場合、LB と UB に空行列を設定してください。
% X(i) に下限がない場合、LB(i) = -Inf と設定し、X(i) に上限がない場合、
% UB(i) = Inf と設定してください。
% 
% X=LSQCURVEFIT(FUN,X0,XDATA,YDATA,LB,UB,OPTIONS) は、オプション 
% (OPTIONS) を設定して実行できます。ここで、OPTIONS は OPTIMSET 関数で
% 設定できる構造体です。詳細は、OPTIMSET を参照してください。Display, 
% TolX, TolFun, DerivativeCheck, Diagnostics, FunValCheck, Jacobian, 
% JacobMult, JacobPattern, LineSearchType, LevenbergMarquardt, MaxFunEvals, 
% MaxIter, DiffMinChange, DiffMaxChange, LargeScale, MaxPCGIter, 
% PrecondBandWidth, TolPCG, OutputFcn, TypicalX  のオプションが使用できます。
% オプション Jacobian を使って、FUN を呼び出し、2 番目の出力引数 J に、
% 点 X でのヤコビ行列を設定できます。[F,J] = feval(FUN,X) の形式で呼び
% 出します。FUN が、X が長さ n の場合、m 要素のベクトルを返す場合、J は、
% m 行 n 列の行列になります。ここで、J(i,j) は、F(i) の x(j) による偏微分
% 係数です (ヤコビ行列 J は、F の勾配を転置したものです)。
%
% [X,RESNORM]=LSQCURVEFIT(FUN,X0,XDATA,YDATA,...)  は、X での残差の 
% 2 ノルム sum {(FUN(X,XDATA)-YDATA).^2} を出力します。
%
% [X,RESNORM,RESIDUAL] = LSQCURVEFIT(FUN,X0,...) は、解 X での残差の値 
% FUN(X,XDATA)-YDATA を出力します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG] = LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) は、
% LSQCURVEFIT の終了状況を示す文字列 EXITFLAG を出力します。EXITFLAG の
% 可能な値と対応する終了状況は、つぎの通りです。
%
%   1  LSQCURVEFIT は、解 X に収束したことを示します。
%   2  X の変化が指定した許容範囲より小さいことを示します。
%   3  目的関数値の変化が指定した許容範囲より小さいことを示します。
%   4  探索方向の大きさが指定した許容誤差より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを示します。
%  -1  最適化が、出力関数で終了していることを示します。
%  -2  境界が矛盾していることを示します。
%  -4  ライン探索は、現在の探索方向に沿って十分に残差を減少できないことを示します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) 
% は、繰り返し回数 OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、
% 使用したアルゴリズム OUTPUT.algorithm、(使用した場合) 共役勾配繰り返し回数 
% OUTPUT.cgiterations 、(使用した場合) 1 次の最適性 OUTPUT.firstorderopt、
% 終了メッセージ OUTPUT.message をもつ構造体 OUTPUT に出力します。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) 
% は、解 X でのラグランジェ乗数 LAMBDA を出力します。LAMBDA.lower に LB を、
% LAMBDA.upper に UB を設定しています。
%
% [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = LSQCURVEFIT(FUN,X0,XDATA,YDATA,...) 
% は、解 X での関数 FUN のヤコビ行列値を出力します。
%
% 例
% FUN は、@ を使って設定することができます。
% 　　   xdata = [5;4;6];          % xdata の例
%        ydata = 3*sin([5;4;6])+6; % ydata の例
%        x = lsqcurvefit(@myfun, [2 7], xdata, ydata)
%
% ここで、myfun は、つぎのように表される MATLAB 関数です。
%
%       function F = myfun(x,xdata)
%       F = x(1)*sin(xdata)+x(2);
%
% FUN は、無名関数としても表現できます。
%
%       x = lsqcurvefit(@(x,xdata) x(1)*sin(xdata)+x(2),[2 7],xdata,ydata)
%
% FUN がパラメータ化された場合、問題に依存したパラメータを指定して無名
% 関数を使用できます。2 番目の引数 c でパラメータ化された関数 myfun を
% 最小化すると仮定します。ここで、mfun はつぎのような M-ファイル関数です。
%
%       function F = myfun(x,xdata,c)
%       F = x(1)*exp(c*xdata)+x(2);
%
% c の指定値に対して最適化するためには、最初に c に設定します。次に、
% x を引数とする無名関数 myfun を定義し、最終的に、LSQCURVEFIT に渡します。
%
%       xdata = [3; 1; 4];           % xdata の例
%       ydata = 6*exp(-1.5*xdata)+3; % ydata の例
%       c = -1.5;                    % パラメータの定義
%       x = lsqcurvefit(@(x,xdata) myfun(x,xdata,c),[5;1],xdata,ydata)
%
%   参考 OPTIMSET, LSQNONLIN, FSOLVE, @, INLINE.


%   Copyright 1990-2006 The MathWorks, Inc.
