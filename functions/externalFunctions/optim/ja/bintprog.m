% BINTPROG  2 値整数計画問題
%
%   BINTPROG は、2 値整数計画問題を解きます。
%   つぎの制約のもとで、f'*X を X に関して最小化します。
%                    A*X <= b,
%                    Aeq*X = beq,
%                    ここで、X の要素は 2 値整数、
%                    すなわち 0 か 1 です。
%
%   X = BINTPROG(f) は、 f'*X を最小化する問題を解きます。ここで、X の
%   要素は、2 値整数です。
%
%   X = BINTPROG(f,A,b) は、線形不等式 A*X <= b のもとで、f'*X を最小化
%   する問題を解きます。ここで、X の要素は、2 値整数です。
%
%   X = BINTPROG(f,A,b,Aeq,beq) は、線形等式 Aeq*X = beq と 線形不等式 
%   A*X <= b のもとで f'*X を最小化する問題を解きます。ここで、X の要素は、
%   2 値整数です。
%
%   X = BINTPROG(f,A,b,Aeq,beq,X0) は、初期値を X0 に設定します。初期値 
%   X0 は、2 値整数で、かつ可解でなければなりません。そうでない場合、
%   無視されます。
%
%   X = BINTPROG(f,A,b,Aeq,beq,X0,OPTIONS) は、オプション (OPTIONS) を
%   設定して実行できます。ここで、OPTIONS は OPTIMSET 関数で設定できる
%   構造体です。詳細は、OPTIMSET を参照してください。BranchStrategy, 
%   Diagnostics, Display, NodeDisplayInterval, MaxIter, MaxNodes, 
%   MaxRLPIter, MaxTime, NodeSearchStrategy, TolFun, TolXInteger, 
%   TolRLPFun のオプションが使用できます。
%
%   [X,FVAL] = BINTPROG(...) は、解 X での目的関数の値 FVAL = f'*X を
%   出力します。
%
%   [X,FVAL,EXITFLAG] = BINTPROG(...) は、BINTPROG の終了状況を示す文字列 
%   EXITFLAG を返します。EXITFLAG の可能な値と対応する終了状況は、つぎの
%   通りです。
%
%      1  BINTPROG は、解 X に収束したことを示します。
%      0  繰り返し回数が最大回数に達していることを示します。
%     -2  可解が見つからなかったことを示します。
%     -4  収束することなく MaxNodes に達したことを示します。
%     -5  収束することなく MaxTime に達したことを示します。
%     -6  LP 緩和問題を解くためのノードで実行された繰り返し回数が、
%         収束することなく MaxRLPIter に達したことを示します。
%
%   [X,FVAL,EXITFLAG,OUTPUT] = BINTPROG(...) は、繰り返し数 OUTPUT.iterations、
%   探索したノードの数 OUTPUT.nodes、実行時間 (秒単位) OUTPUT.time、
%   使用したアルゴリズム OUTPUT.algorithm、分岐法 OUTPUT.branchStrategy、
%   ノード検索方法 OUTPUT.nodeSrchStrategy、終了メッセージ OUTPUT.message 
%   をもつ構造体 OUTPUT を出力します。
%
%   例
%     f = [-9; -5; -6; -4];
%     A = [6 3 5 2; 0 0 1 1; -1 0 1 0; 0 -1 0 1];
%     b = [9; 1; 0; 0];
%     X = bintprog(f,A,b)


%   Copyright 1990-2006 The MathWorks, Inc.
