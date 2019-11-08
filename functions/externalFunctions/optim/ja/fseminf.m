% FSEMINF   半無限制約付き最適化問題
%
% FSEMINF は、つぎの形式の問題を解きます。
%
%         { F(x) | C(x) <= 0 , Ceq(X) = 0 , PHI(x,w) <= 0 } 
% 区間内に存在するすべての w に対して、F(x) を x に関して最小化します。
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON) は、X0 を初期値として、関数 SEMINFCON 
% で規定された NTHETA 個の半無限制約のもとで、関数 FUN 最小化する X を
% 求めます。関数 FUN は、X をベクトル入力として、X で計算されるスカラの
% 関数値 F を出力します。関数 SEMINFCON は、X, S を非線形入力として、
% 非線形不等式制約ベクトル C、非線形等式制約ベクトル Ceq、および、ある
% 区間全体で評価される NTHETA 個の半無限不等式制約行列 PHI_1, PHI_2, ..., 
% PHI_NTHETA を出力します。S は、推奨するサンプル間隔で、使用してもしなくても
% 構いません。
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B) は、線形不等式制約 A*X <=  B も
% 満足させようとします。
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq) は、線形等式制約 
% Aeq*X = Beq も満足させようとします (不等式が存在しない場合、A = [] と 
% B = [] を設定してください)。
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq,LB,UB) は、設計変数 X 
% の上下限の範囲を与えることが可能です。この場合、LB <= X <= UB の範囲内で
% 最適解が探索されます。範囲の制約がない場合、LB と UB に空行列を設定
% してください。X(i) に下限がない場合、LB(i) = -Inf と設定し、X(i) に
% 上限がない場合、UB(i) = Inf と設定してください。
%
% X = FSEMINF(FUN,X0,NTHETA,SEMINFCON,A,B,Aeq,Beq,LB,UB,OPTIONS) は、
% オプション (OPTIONS) を設定して実行できます。ここで、OPTIONS は OPTIMSET 
% 関数で設定できる構造体です。詳細は、OPTIMSET を参照してください。Display, 
% TolX, TolFun, TolCon, DerivativeCheck, Diagnostics, FunValCheck, GradObj, 
% MaxFunEvals, MaxIter, DiffMinChange, DiffMaxChange, OutputFcn, TypicalX 
% のオプションを使用できます。オプション GradObj を使って、FUN を呼び出し、
% 2 番目の出力引数 G に点 X での偏微分係数 df/dX を設定できます。
% [F,G] = feval(FUN,X) の形式で呼び出します。
%
% [X,FVAL] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) は、解 X での目的関数
% FUN の値を出力します。
%
% [X,FVAL,EXITFLAG] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) は、FSEMINF の
% 終了状況を示す文字列 EXITFLAG を出力します。EXITFLAG の可能な値と対応する
% 終了状況は、つぎの通りです。
% 
%   1  FSEMINF は解 X に収束したことを示します。
%   4  探索方向の大きさが指定した許容誤差より小さく、制約違反が options.TolCon 
%      より小さいことを示します。
%   5  探索方向に沿った関数の勾配の大きさが許容範囲より小さく、制約違反が 
%      options.TolCon より小さいことを示します。
%   0  関数計算の回数、あるいは繰り返し回数が最大回数に達していることを示します。
%  -1  最適化が、出力関数で終了していることを示します。
%  -2  可解が見つからなかったことを示します。
% 
% [X,FVAL,EXITFLAG,OUTPUT] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...) は、
% 繰り返し数 OUTPUT.iterations、関数の評価回数 OUTPUT.funcCount、
% 最終ステップのノルム OUTPUT.stepsize、最終のライン探索のステップ長 
% OUTPUT.lssteplength、使用したアルゴリズム OUTPUT.algorithm、1 次の最適性 
% OUTPUT.firstorderopt、終了メッセージ  OUTPUT.message をもつ構造体 
% OUTPUT を出力します。
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = FSEMINF(FUN,X0,NTHETA,SEMINFCON,...)
% は、解でのラグランジェ乗数を出力します。LAMBDA 構造体は、LAMBDA.lower
% に LB を、LAMBDA.upper に UB を、LAMBDA.ineqlin に線形不等式を、LAMBDA.eqlin 
% に線形等式を、LAMBDA.ineqnonlin に非線形不等式を、LAMBDA.eqnonlin に
% 非線形等式を設定します。
% 
% 例
% FUN と SEMINFCOM は、@ を使って、設定することができます。
%
%        x = fseminf(@myfun,[2 3 4],3,@myseminfcon)
%
% ここで、myfun は、つぎのように表わされる MATLAB 関数です。
% 
%    function F = myfun(x)
%    F = x(1)*cos(x(2))+x(3)^3:
% 
% また、myseminfcon は、つぎのように表わされる MATLAB 関数です。
% 
%       function [C,Ceq,PHI1,PHI2,PHI3,S] = myseminfcon(X,S)
%       C = [];     % C と Ceq を計算するコード:
%                   % 制約がない場合は空行列にします。
%       Ceq = [];
%       if isnan(S(1,1))
%          S = [...] ; % S は ntheta行×2列
%       end
%       PHI1 = ... ;       % PHI を計算するコード
%       PHI2 = ... ;
%       PHI3 = ... ;
% 　　　　　
%   参考 OPTIMSET, @, FGOALATTAIN, LSQNONLIN.


%   Copyright 1990-2006 The MathWorks, Inc.
