************
インストール
************

動作環境
========

ProteinDFはいくつかの形式で配布されています。
READMEファイルが添付されている場合は、まずREADMEファイルをご覧下さい。


ProteinDFの動作には以下の動作環境が必要です。

* UNIX/Linux オペレーティングシステム
* Cランタイムライブラリ(通常はOSに同梱されています)
* bash
* BLAS, LAPACKライブラリ
* MPI実行環境(並列版のみ)
* 分散行列演算ライブラリ(ScaLAPACK)


ハードウェアとオペレーティングシステム
--------------------------------------

POSIX準拠のコンピュータシステムで動作します。
現在、以下の計算機システムで動作確認が行われています。

* SGI社製 Altix 3000シリーズ
* Cray社製 XT-5, XT-6

一般的なx86 PC Linuxでも動作します。


メモリとディスク
----------------

計算モデルの大きさに応じて必要なメモリ量が異なります。
また、並列計算を行う場合、
行列演算にLAPACKを使用する場合とScaLAPACKを使用する場合でも、
1ノードあたりに必要なメモリ量が変わります。
LAPACKを使用する場合は、各ノードに搭載されているメモリ容量が計算可能なサイズの上限です。
一方、ScaLAPACKを使用する場合は、全ノードで計算領域を分散保持しますので、
全ノードのメモリ容量が計算可能サイズの上限になります。
ただし、この他にも計算可能サイズを決定する要因がありますので
目安としてください。

.. warning::

  32bit OSでは扱えるメモリサイズ、ファイルサイズなど幾つかの制限事項がある場合があります。


pythonモジュール
----------------

ProteinDFの動作には、
いくつかのpythonスクリプトが用意されています。
これらpythonスクリプトは、ProteinDFの実行そのものには必要ありませんが、
計算結果の解析用として用意されています。
これらpythonスクリプトの動作には以下のソフトウェア(モジュール)が必要です。
これらのソフトウェアの環境構築は、それぞれのシステムの方針に従って下さい。

* python(version 2.5以上)
* argparseモジュール
* numpyモジュール
* matplotlibモジュール
* MessagePackモジュール
* YAMLモジュール



インストールと準備
==================

配布パッケージによって、インストールの形態が異なります。
READMEファイルが添付されている場合は、READMEファイルの指示に従って下さい。


.. index:: 環境変数

環境変数
--------


ProteinDFの実行には、以下の環境変数を設定する必要があります。
利用する環境に応じて、適切に環境変数を設定してください。

.. index:: PDF_HOME

PDF_HOME
^^^^^^^^

ProteinDFパッケージをコピーしたディレクトリを指定します。


例えば、/usr/local/ProteinDF に本パッケージをコピーしたとき、
ログインシェルにbashを利用している場合は.bashrcに以下を追加してください。

.. code-block:: bash

  export PDF_HOME=/usr/local/ProteinDF


.. index:: OMP_NUM_THREADS

OMP_NUM_THREADS
^^^^^^^^^^^^^^^

OpenMPにおける最大スレッド数を設定します。
ビルド時にOpenMPを有効にする必要があります。


.. index:: OMP_SCHEDULE

OMP_SCHEDULE
^^^^^^^^^^^^

OpenMPにおける並列スケジュールのタイプとチャンクサイズを設定できます。
ビルド時にOpenMPを有効にする必要があります。



ソースからのビルド
------------------

.. index: configure

cmakeの実行
^^^^^^^^^^^^^^^

ProteinDFではMakefileの作成にcmakeを利用しています。
任意のディレクトリにて、ソースディレクトリを指定してcmakeを実行します。
例えばtarballを展開したディレクトリの直下にbuildディレクトリを作成する場合、
以下のように実行します。

.. code-block:: bash

  $ mkdir build
  $ cd build
  $ cmake ..


.. note::

  cmakeは自動的にビルド環境を調査し、ライブラリの場所を設定します。
  調査結果はcmakeの実行時に出力されます。
  出力結果をファイルに保存して後に参照したい場合は、teeコマンドを利用して
  ``$ cmake .. 2>&1 | tee out.cmake``
  のように実行します。


以下によく用いられる変数を示します。
詳しくは ``cmake -L`` または ``cmake -LA`` をご覧ください。


* ``--prefix=location``

プログラムのインストール先を指定します。
デフォルトは/usr/localです。
ユーザーのホームディレクトリなどにインストールするときに用いられます。


* ``--enable-parallel``

逐次版に加えて並列版プログラムも作成します。


* ``--with-blas=location``

BLASライブラリの場所を指定します。


* ``--with-lapack=location``

LAPACKライブラリの場所を指定します。


* ``--with-scalapack=location``

ScaLAPACKライブラリの場所を指定します。


* 環境変数CC,CXX,MPICXX

ビルドに用いるC/C++コンパイラを指定します。
MPIライブラリをリンクする場合は、
mpicxxなど計算機システムに応じたコンパイラを使用してください。


* 環境変数CFLAGS,CXXFLAGS

C/C++コンパイラに渡すオプションを指定します。
OpenMPなどのコンパイラへの指示は、この環境変数に指定してください。

* 環境変数LIBS

その他ビルドに必要なライブラリを指定します。

* 環境変数BLAS_LIBS

BLASライブラリを指定します。

* 環境変数LAPACK_LIBS

LAPACKライブラリを指定します。

* 環境変数SCALAPACK_LIBS

ScaLAPACKライブラリを指定します。


makeの実行
^^^^^^^^^^

configureを実行し、Makefileが作成されたことを確認してください。
Makefileが作成されたならば、makeを実行します。


.. code-block:: bash

  $ make 2>&1 | tee out.make


インストールの実行
^^^^^^^^^^^^^^^^^^

makeの実行した後、実行ファイル・データを所定のパスにインストールします。


.. code-block:: bash

  $ make install 2>&1 | tee out.make_install


インストールが実行されると、以下のファイルがコピーされます。


* ${PDF_HOME}/bin/PDF.x
* ${PDF_HOME}/bin/PPDF.x
* ${PDF_HOME}/data/basis2


うまく行かない場合
^^^^^^^^^^^^^^^^^^

環境によってはスクリプトが実行できない、ビルドできないなどの問題が発生するかもしれません。
その場合は、各スクリプトの出力(上記の操作の場合、out.configure, out.make, out.make_install)をよくチェックしてください。


