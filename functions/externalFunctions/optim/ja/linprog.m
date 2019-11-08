% LINPROG  線形計画問題
%
% X = LINPROG(f,A,b) は、つぎの形式の線形計画問題を解きます。
% A*x <= b のもとで、f'*x を x に関して最小化します。
%  
% X = LINPROG(f,A,b,Aeq,beq) は、等式制約 Aeq*X = Beq を考慮して、上記の
% 問題を解きます。
% 
% X = LINPROG(f,A,b,Aeq,beq,LB,UB) は、設計変数 X の上下限の範囲を与える
% ことが可能です。この場合、LB <= X <= UB の範囲内で最適解が探索されます。
% 範囲の制約がない場合、LB と UB に空行列を設定してください。X(i) に
% 下限がない場合、LB(i) = -Inf と設定し、X(i) に上限がない場合、UB(i) = Inf 
% と設定してください。
%
% X = LINPROG(f,A,b,Aeq,beq,LB,UB,X0) は、初期値を X0 とします。この
% オプションは、active-set アルゴリズムを使用するときのみ利用できます。
% デフォルトの内点アルゴリズムは、全ての空でない初期値を無視します。
%
% X=LINPROG(f,A,b,Aeq,beq,LB,UB,X0,OPTIONS) は、オプション (OPTIONS) 
% を設定して実行できます。ここで OPTIONS は OPTIMSET 関数で設定できる
% 構造体です。詳細は、OPTIMSET を参照してください。Display, Diagnostics, 
% TolFun, LargeScale, MaxIter のオプションが使用できます。LargeScale が 
% 'off' の場合、'final' と 'off' のみが、Display に対して使われます 
% (LargeScale が 'on' の場合、'iter' を使用できます)。
%
% [X,FVAL] = LINPROG(f,A,b) は、解 X での目的関数の値 FVAL = f'*X を
% 出力します。
%
% [X,FVAL,EXITFLAG] = LINPROG(f,A,b) は、LINPROG の終了状況を示す文字列 
% EXITFLAG を出力します。EXITFLAG の可能な値と対応する終了状況は、つぎの
% 通りです。
%
%   1  LINPROG は、解 X に収束したことを示します。
%   0  繰り返し回数が最大回数に達していることを示します。
%  -2  可解が見つからなかったことを示します。
%  -3  問題に制約がないことを示します。
%  -4  アルゴリズムの評価中に NaN の値があったことを示します。
%  -5  双対、主問題の両方が不可解であることを示します。
%  -7  探索方向の大きさが小さくなりすぎていることを示します。これ以上
%      計算を行なうことができません。
%
% [X,FVAL,EXITFLAG,OUTPUT] = LINPROG(f,A,b) は、繰り返し回数 OUTPUT.iterations、
% 使用したアルゴリズム OUTPUT.algorithm、(使用した場合) 共役勾配繰り返し回数 
% OUTPUT.cgiterations、終了メッセージ OUTPUT.message をもつ構造体 OUTPUT 
% を出力します。
%
% [X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = LINPROG(f,A,b) は、解 X でのラグランジェ
% 乗数 LAMBDA を出力します。LAMBDA 構造体は、LAMBDA.lower に LB を、
% LAMBDA.upper に UB を、LAMBDA.ineqlin に線形不等式を、LAMBDA.eqlin に
% 線形等式を、LAMBDA.ineqnonlin に非線形不等式を、LAMBDA.eqnonlin に
% 非線形等式を設定します。
%   
% 注意 ： LINPROG の大規模問題 (デフォルト) は、主双対法を使用します。
% 主問題と双対問題は、共に収束に対して不可解でなければなりません。主問題、
% 双対問題、あるいは両方のいずれかが可解であるというメッセージは、適切に
% 与えられます。主問題の標準的な形式は、つぎの通りです。
% 
%      A*x = b, x >= 0 のもとで、f'*x を最小化します。
% 
% 双対問題は、つぎの形式になります。
% 
%      A'*y + s = f, s >= 0 のもとで、b'*y を最大化します。


%   Copyright 1990-2006 The MathWorks, Inc.
