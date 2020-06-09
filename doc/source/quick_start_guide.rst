**********************
クイックスタートガイド
**********************

ここではアミノ酸(グリシン)の一点計算を1プロセスで計算します。

.. index:: fl_Input
.. index:: fl_Table
.. index:: fl_Work

計算用ディレクトリの用意
========================

ProteinDFでは行列・ベクトルなどのサイズ・数ともに大きな中間ファイルを
ディスクに書き込みます。
ProteinDFの実行に際し、あらかじめ計算に使用するディレクトリを作成します。
作成すべきディレクトリは

* fl_Work

です。

.. note::

  ``${PDF_HOME}/bin/pdf-setup`` を実行すると、計算に必要な環境を整えます。


.. warning::

  ディスクの書き込みに失敗すると、プログラムが異常終了することがあります。
  特にMPI並列計算を行う場合は、すべてのノードから書き込みができることを確認してください。


入力ファイルの準備
==================

.. index:: fl_Userinput

以下の内容のテキストファイルを作成し、
``fl_Userinput``
というファイル名で保存します。


.. code-block:: text
  
  >>>>MAIN
        step-control    = [create integral guess scf]
        cut-value       = 1.0e-10
        scf-start-guess = harris
        max-iteration   = 100
        method  = rks
        method/rks/electron-number              = 40
        method/rks/occlevel                     = [ 1 - 20 ]
        orbital-independence-threshold          = 0.007
        convergence/type                        = density
        convergence/threshold                   = 1e-4
        convergence/threshold-energy            = 1e-5
        scf-acceleration                        = damping
        scf-acceleration/damping/damping-factor = 0.65 
        xc-potential                            = b3lyp
        scf-acceleration/damping/damping-type   = density_matrix
  
  >>>>MOLECULE
        geometry/cartesian/unit = angstrom
        geometry/cartesian/input        = {
                N        -1.888000        0.035000       -0.211000
                H        -1.766000        0.945000        0.189000
                H        -1.817000        0.099000       -1.205000
                C        -0.758000       -0.730000        0.287000
                H        -0.893000       -0.915000        1.372000
                H        -0.720000       -1.725000       -0.200000
                C         0.529000        0.065000        0.064000
                O         0.520000        1.294000        0.114000
                O         1.742000       -0.451000       -0.186000
                H         1.692000       -1.400000       -0.203000
        }end
  
        basis-set/orbital       = {
                H = "O-HYDROGEN (41) DZVP"
                O = "O-OXYGEN (621/41) by FS"
                C = "O-CARBON (621/41) by FS"
                N = "O-NITROGEN (621/41) by FS"
        }end
  
        basis-set/density-auxiliary     = {
                H = "A-HYDROGEN (4,1;4,1) from deMon"
                O = "A-OXYGEN (7/2;7/2) by FS"
                C = "A-CARBON (7/2;7/2) by FS"
                N = "A-NITROGEN (7/2;7/2) by FS"
        }end
  
        basis-set/exchange-auxiliary    = {
                H = "A-HYDROGEN (4,1;4,1) from deMon"
                O = "A-OXYGEN (7/2;7/2) by FS"
                C = "A-CARBON (7/2;7/2) by FS"
                N = "A-NITROGEN (7/2;7/2) by FS"
        }end


.. note::

  ProteinDFの起動オプションにより、
  入力ファイルを変更することができます。


ProteinDFの実行
===============

環境変数PDF_HOMEを適切に設定した後、
ProteinDF(逐次版)を実行します。

.. code-block:: bash

  % ${PDF_HOME}/bin/PDF.x


正常に終了した場合は、コマンドプロンプトに戻ります。


結果の表示
==========

計算結果はファイルに出力されます。
以下に例を示します。


.. note::

  出力ファイルの場所は、ProteinDFの起動オプションにより変更できます。


はじめにProteinDFのバージョン、ならびに並列数(MPIプロセス数、OpenMPスレッド数)が
表示されます。
意図した通りに実行されているか確認してください。

.. code-block:: none
   
   [0:2012/**/07 17:17:02:INFO] **************************************
   [0:2012/**/07 17:17:02:INFO] ProteinDF version 20xx.x:xxxx (serial)
   [0:2012/**/07 17:17:02:INFO] 
   [0:2012/**/07 17:17:02:INFO]  OpenMP threads: 12
   [0:2012/**/07 17:17:02:INFO] 

``step_control`` に記載されている手順に従い、計算が実行されます。
ログの左側に出力日時が記載されます。

.. code-block:: none
   
   ===============================================
    >>>> INTEGRAL
   ===============================================
    >>>> Hpq
   
    ...
   
   ===============================================
    >>>> GUESS
   ===============================================
   
    ...
   
   ===============================================
    >>>> SCF
   ===============================================


エネルギー情報は ``Total Energy`` ブロックに出力されます。

.. code-block:: none
   
   ------------------------------------------------
    >>>> Total Energy
   ------------------------------------------------
    Ts+Vn          =        -745.7264230071891689
    E_J[Rho, Rho~] =         629.4907739256434525
    E_J[Rho~,Rho~] =        -311.1519936178789294
    E_xc(pure)     =         -27.8222685895087842
    E_K            =          -6.8830892822011522
    E_nuclei       =         179.8930288412503558
    TE             =        -282.1999717298841688
   ------------------------------------------------

その他、計算サイズ、カットオフ情報が随時出力されます。

計算が正常に終了すると、CPU時間と経過時間が出力されます。

.. code-block:: none
   
   ************************************************
    ProteinDF successful completion
    CPU_TIME:         3454 sec
    ELAPS_TIME:        542 sec
   ************************************************
