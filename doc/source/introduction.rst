.. -*- coding: utf-8; -*-

****
概要
****

ProteinDFはタンパク質の精密な全電子カノニカル分子軌道計算の実行・解析に焦点をあてたプログラムです。



特徴
====

ProteinDFには以下の特徴があります。

* Gauss型基底関数を用いた密度汎関数法プログラム

  * Hartree-Fock計算 および 純粋密度汎関数計算、混合密度汎関数計算

* オブジェクト指向言語C++で開発

* タンパク質をはじめとする巨大分子の基底状態全電子カノニカル分子軌道計算

  * RI(Resolution of Identity)法による分子積分の高速化
  * コレスキー分解法による精密化と高速化の両立
  * 大規模行列の分散保持による大規模計算
  * MPI / OpenMPによるハイブリッド並列化による並列計算

* 巨大分子の解析

  * 各原子座標におけるエネルギー勾配の解析的計算
  * Mulliken電子密度解析
  * 分子軌道・電子密度・静電ポテンシャルの計算


計算手法
========

タンパク質の構造を丸ごと取り扱い、機能を定量的に予測できる計算方法は、
電子相関を効果的に取り込むことができるKohn-Sham方程式を解く密度汎関数法が適しています。
ProteinDFは基底関数にガウス型関数を使用するKohn-Sham-Roothaan法に基づく密度汎関数計算プログラムで、
化学では標準的な分子軌道計算法です。
密度汎関数法の解説や、大規模計算技術の詳細は他書に譲り、
ここでは簡単にProteinDFの概要を述べます。

大規模なカノニカル計算を達成するため、
計算高速化法としてRI (Resolution of Identity)法やコレスキー分解法を採用しています。
計算方法は制限付き、非制限付き、および制限付き開殻Kohn-Sham(Restricted, Unrestricted and Restricted Open-shell KS) 法をサポートしています。
デフォルトの交換相関ポテンシャルはハイブリッド汎関数B3LYPで、
局所汎関数SVWN、ならびに一般化密度勾配補正のBLYPなどが用意されています。
なお、SVWNに関しては、交換相関ポテンシャル計算にもRI法を用いて計算する方法も選択することができます。
これらの方法で、エネルギー勾配計算が用意されていますが、
これを用いた構造最適化機能や量子分子動力学法計算、
物理量の計算には、これらを実行するための外部プログラムと連携を行う必要があります。
基底関数と補助基底関数はガウス型基底関数を用いる密度汎関数法に最適化したものを
basis2 ファイルに用意しています。追加も可能です。
ProteinDF のデフォルトセットはValence Double(VD) 相当のものを指定してあります。
詳しくはbasis2 ファイルをご覧下さい。
ProteinDF は他の分子軌道法プログラムと同様にキーワード方式で実行が指定されます。
上記の説明中にも現れましたが、ほぼすべてのキーワードはデフォルト値が設定されています。
詳しくは付録を参考にして下さい。


計算サイズ
==========

ProteinDF の最大の特徴はその計算サイズです。
数十残基からなるタンパク質の全電子計算が現在のPCで無理なく実行できる理由は、
私たちが行ってきたさまざまな工夫に起因するものと自負していますが、
もっとも大きな要因はメモリ仕様に関するプログラミングの取り決めにあります。
ProteinDF では動的なメモリ管理により、
基底関数の次元をもつ行列を数個分のみを使用して全実行が行われます。
この事実から逆に、使用するマシンのメモリサイズから最大計算サイズが算出できます。
ちなみに、基底関数の数をN :sub:`orb` とすると、行列1つ分のメモリは8×N :sub:`orb` × N :sub:`orb` バイトと概算できます。
ただし、通常OS や転送などが使用するメモリ領域がありますし、
特に動的なメモリ消去において(たとえプログラムに明記していても)、OS が勝手なタイミングで行う場合があります。
また、メモリを最大限に利用するため、中間出力としてディスク媒体を使用しています。
計算中は計算サイズの行列を最低200個分保持できるよう、ディスクの空きを確保して下さい。
ProteinDF では安全のため中間出力をすべて残しています。


その他
======

ProteinDF はペプチド鎖からなるタンパク質の全電子計算に適化しています。
一方、ヘテロ分子を含む計算など、より複雑な計算には自動計算法もGUI の対応もまだ不十分な状態です。
このような計算を行うには、本手法やプログラムを深く理解し、
場合によっては手作業が必要となると思います。
詳しくはリファレンスマニュアルを参考にして下さい。
最後に、これまで私達がペプチド鎖からなるタンパク質の全電子計算を通して得た経験のうち重要と思われるものを以下に列挙します。


タンパク質サイズの経験則
------------------------

タンパク質の密度汎関数計算とはどの程度の規模の計算なのかを確認することは、
初めに行うべき最も大切なことの一つです。
Roothaan法ではKohn-Sham方程式の固有値を求める問題を行列方程式の固有値解法に置き換えます。
そのため、行列のサイズが計算の規模を表す指標になります。
昔から、タンパク質ではアミノ酸残基の総数と分子量との間に次の経験則が知られています。
 
(タンパク質の分子量) = 110 × (アミノ酸残基数)

この関係式は、タンパク質がペプチド鎖でできていること、
アミノ酸をつくる水素、炭素、窒素、酸素、硫黄といった原子の組成比率が巨大な分子で平均化されるためによく合います。
この経験則を発展させて、タンパク質のアミノ酸残基数(N\ :sub:`res`)、原子数(N\ :sub:`atom`)、
電子数(N\ :sub:`ele`)、軌道数(N\ :sub:`orb`)(すなわち行列の次元)との間の比例関係を導くことができます。
例えばProteinDF のデフォルトセット(VD 相当) ではおおよそ以下の関係が成り立ちます。

N\ :sub:`res` : N\ :sub:`atom` : N\ :sub:`ele` : N\ :sub:`orb` = 1 : 20 : 70 : 100

ここで、N\ :sub:`atom` のうち半分は水素原子であることも経験的にわかっています。
これは計算サイズを知るうえで有用な関係式です。
もちろん、N\ :sub:`orb` は使用する基底関数等のセットに依存しますが、
異なるセットでも同様の比例関係式を簡単に見積もることができます。
計算を行う前に、使用するセットにおける比例関係式を導いて、安全な計算サイズを把握して下さい。


タンパク質構造の歪み
--------------------

タンパク質の立体構造はProtein Data Bank(PDB) から入手するのが一般的です。
PDBは、X 線構造解析や中性子散乱、多次元NMR などで決定したタンパク質構成原子の3 次元座標データ群です。
しかし、PDB の座標データには実験そのものの性格上、あるいはその後のデータチューニングの性質上、
異常な構造の歪みを持っているものが多数あります。
これは密度汎関数法計算において特に悪影響を与えます。
また、現在の計算機資源では、全電子カノニカル分子軌道法による密度汎関数法そのもので構造最適化することはまだ現実的ではありません。
計算の際には分子力学法、MOZYME などの半経験法、QM/MMやONIOM法、
フラグメント分子軌道法などで、最適化されることをお勧めします。
また、タンパク質の構造歪みのチェックには例えばPROCHECK などを参考にして下さい。


タンパク質表面の取り扱い
------------------------

化学分子では孤立系のシミュレーションが標準的に行われていますが、
これまでの経験では、密度汎関数法では表面に解離基が多数存在する水溶性タンパク質を真空中で解くことはできません。
これは密度汎関数法の精度の良さを示す1 つの事実であると考えられますが、
一方で計算者側からすれば大変厄介な事実でもあります。
タンパク質の表面がどのようになっているのかは、
今でも最先端の研究課題です。
最も良い方法はタンパク質の周りに溶媒分子(水分子やバッファーイオン)を適切に相当数配置させることですが、
全て量子的に取り扱うと計算サイズが飛躍的に膨れ上がってしまいます。
そのような計算はチャレンジングで大変興味がそそられますが、
通常の計算では、タンパク質周りドロップレット状に古典的な水分子やカウンターイオンを配置させる、
といった便法を使用するなどの対処が必要です。


