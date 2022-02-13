* pdf_param_path

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "pdf_param_path"
  "explanation", ""
  "default", "pdfparam.mpac"
  "values", ".*"


* step_control

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "step_control"
  "explanation", ""
  "default", "create"
  "values", "((create|integral|guess|scf|force)\s*)+"


* comment

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "comment"
  "explanation", ""
  "default", ""
  "values", ".*"


* linear_algebra_package

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "linear_algebra_package"
  "explanation", ""
  "default", "lapack"
  "values", "(lapack|scalapack)"


* scalapack_block_size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scalapack_block_size"
  "explanation", ""
  "default", "64"
  "values", "(\d+)"


* save_distributed_matrix_to_local_disk

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "save_distributed_matrix_to_local_disk"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* local_temp_path

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "local_temp_path"
  "explanation", ""
  "default", ""
  "values", ".*"


* memory_size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "memory_size"
  "explanation", ""
  "default", "1GB"
  "values", "\d+\s*(MB|GB)"


* use_mmap

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "use_mmap"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* mapfile_size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "mapfile_size"
  "explanation", ""
  "default", "auto"
  "values", "(\d+\s*(MB|GB)|auto)"


* mapfile_basename

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "mapfile_basename"
  "explanation", ""
  "default", "pdf"
  "values", ".*"


* work_on_disk

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "work_on_disk"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* parallel_processing_type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "parallel_processing_type"
  "explanation", ""
  "default", "divide_and_conquer"
  "values", "(divide_and_conquer|DC|master_slave|MS)"


* cleanup

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "cleanup"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* show_keyword

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "show_keyword"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* show_input

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "show_input"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* show_coordinates

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "show_coordinates"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* show_orbital_basis

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "show_orbital_basis"
  "explanation", ""
  "default", "gamess"
  "values", "(gamess|amoss|none)"


* guess

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess"
  "explanation", ""
  "default", "rho"
  "values", "(rho|file_rho|lcao|density|core|huckel|harris)"


* guess/normalize_density_matrix

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/normalize_density_matrix"
  "explanation", ""
  "default", "true"
  "values", "string"


* cut_value

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "cut_value"
  "explanation", ""
  "default", "1.0E-10"
  "values", "(double > 0)"


* block-size

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "block-size"
  "explanation", ""
  "default", "1024000"
  "values", "(integer > 0)"


* charge-extrapolate-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "charge-extrapolate-number"
  "explanation", ""
  "default", "1"
  "values", "(integer >= 1)"


* orbital-overlap-correspondence

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-overlap-correspondence"
  "explanation", ""
  "default", "off"
  "values", "on/off"


* orbital-overlap-correspondence-first

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-overlap-correspondence-first"
  "explanation", ""
  "default", "off"
  "values", "on/off"


* orbital-overlap-correspondence-method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-overlap-correspondence-method"
  "explanation", ""
  "default", "mo-overlap"
  "values", "mo-overlap/mo-projection"


* orbital-overlap-correspondence-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-overlap-correspondence-number"
  "explanation", ""
  "default", "3"
  "values", "integer"


* summary

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "summary"
  "explanation", ""
  "default", "none"
  "values", "(none|convergence|every-scf)"


* analyze_population

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "analyze_population"
  "explanation", ""
  "default", "none"
  "values", "(none|convergence|every-scf)"


* max_iteration

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "max_iteration"
  "explanation", ""
  "default", "100"
  "values", "(integer)"


* method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method"
  "explanation", ""
  "default", "rks"
  "values", "(rks|uks|roks)"


* method/rks/occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/rks/occlevel"
  "explanation", ""
  "default", ""
  "values", "(array of integer >= 0)"


* method/rks/electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/rks/electrons"
  "explanation", ""
  "default", ""
  "values", "(integer >= 2)"


* method/uks/alpha_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/uks/alpha_occlevel"
  "explanation", ""
  "default", ""
  "values", "(array of integer >= 0)"


* method/uks/alpha_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/uks/alpha_electrons"
  "explanation", ""
  "default", ""
  "values", "(integer >= 1)"


* method/uks/beta_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/uks/beta_occlevel"
  "explanation", ""
  "default", ""
  "values", "(array of integer >= 0)"


* method/uks/beta_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/uks/beta_electrons"
  "explanation", ""
  "default", ""
  "values", "(integer >= 1)"


* method/roks/closed_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/roks/closed_occlevel"
  "explanation", ""
  "default", ""
  "values", "(array of integer >= 0)"


* method/roks/closed_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/roks/closed_electrons"
  "explanation", ""
  "default", ""
  "values", "(integer >= 1)"


* method/roks/open_occlevel

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/roks/open_occlevel"
  "explanation", ""
  "default", ""
  "values", "(array of integer >= 1)"


* method/roks/open_electrons

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "method/roks/open_electrons"
  "explanation", ""
  "default", ""
  "values", "(integer >= 0)"


* save_diff_density_matrix

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "save_diff_density_matrix"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* use_matrix_cache

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "use_matrix_cache"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* force_loading_from_disk

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "force_loading_from_disk"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* show_cache_report

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "show_cache_report"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* disk-utilization

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "disk-utilization"
  "explanation", ""
  "default", "no"
  "values", "yes/no"


* update_method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "update_method"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* orbital-independence-threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-independence-threshold"
  "explanation", ""
  "default", "0.007"
  "values", "(real >= 0)"


* orbital-independence-threshold/canonical

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-independence-threshold/canonical"
  "explanation", ""
  "default", ""
  "values", "(real >= 0)"


* orbital-independence-threshold/lowdin

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "orbital-independence-threshold/lowdin"
  "explanation", ""
  "default", ""
  "values", "(real >= 0)"


* convergence/type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "convergence/type"
  "explanation", ""
  "default", "density_matrix"
  "values", "(fock|density|dcoef)"


* convergence/threshold/rms

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "convergence/threshold/rms"
  "explanation", ""
  "default", "0.001"
  "values", "(real > 0)"


* convergence/threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "convergence/threshold"
  "explanation", ""
  "default", "1e-3"
  "values", "(real > 0)"


* convergence/threshold_energy

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "convergence/threshold_energy"
  "explanation", ""
  "default", "1e-4"
  "values", "(real > 0)"


* scf_acceleration

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration"
  "explanation", ""
  "default", "damping"
  "values", "(damping|anderson|diis)"


* scf_acceleration/damping/damping_factor

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/damping/damping_factor"
  "explanation", ""
  "default", "0.85"
  "values", "(real)"


* scf_acceleration/daming/start_number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/daming/start_number"
  "explanation", ""
  "default", "0"
  "values", "(int)"


* scf_acceleration/damping/damping_type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/damping/damping_type"
  "explanation", ""
  "default", "density_matrix"
  "values", "(fock|density_matrix|density|dcoef)"


* scf_acceleration/oda/start

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/oda/start"
  "explanation", ""
  "default", "10"
  "values", "(int)"


* scf_acceleration/anderson/start_number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/anderson/start_number"
  "explanation", ""
  "default", "3"
  "values", "(int)"


* scf_acceleration/anderson/damping_factor

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/anderson/damping_factor"
  "explanation", ""
  "default", "0.50"
  "values", "(real)"


* scf_acceleration/diis/start

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/diis/start"
  "explanation", ""
  "default", "3"
  "values", "(integer >= 0)"


* scf_acceleration/diis/last_items

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf_acceleration/diis/last_items"
  "explanation", ""
  "default", "3"
  "values", "(integer)"


* level_shift

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "level_shift"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* level_shift/start_iteration

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "level_shift/start_iteration"
  "explanation", ""
  "default", "1"
  "values", "(int)"


* level_shift/closed_mo

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "level_shift/closed_mo"
  "explanation", ""
  "default", "0.00"
  "values", "(real)"


* level_shift/open_mo

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "level_shift/open_mo"
  "explanation", ""
  "default", "0.00"
  "values", "(real)"


* level_shift/virtual_mo

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "level_shift/virtual_mo"
  "explanation", ""
  "default", "0.00"
  "values", "(real)"


* xc_functional

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "xc_functional"
  "explanation", ""
  "default", "svwn"
  "values", "(svwn~|svwn|blyp|b3lyp)"


* grid_free

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid_free"
  "explanation", ""
  "default", "false"
  "values", "(yes|no)"


* xc-potential/grid-type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "xc-potential/grid-type"
  "explanation", ""
  "default", "sg-1"
  "values", "fine/medium-fine/medium/coarse/sg-1"


* xc-update

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "xc-update"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* xc-density-threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "xc-density-threshold"
  "explanation", ""
  "default", "1.0E-16"
  "values", "real"


* grid/weight_threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/weight_threshold"
  "explanation", ""
  "default", "1.0E-16"
  "values", "real"


* grid/num_of_radial_shells

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/num_of_radial_shells"
  "explanation", ""
  "default", "75"
  "values", "int"


* grid/num_of_angular_points

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/num_of_angular_points"
  "explanation", ""
  "default", "302"
  "values", "int"


* grid/radial_quadorature_method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/radial_quadorature_method"
  "explanation", ""
  "default", "GC"
  "values", "(GC|Gauss-Chebyshev|EM|Eular-Maclaurin)"


* grid/GC_mapping_type

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/GC_mapping_type"
  "explanation", ""
  "default", "TA"
  "values", "(Becke|TA|KK)"


* grid/partitioning_method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/partitioning_method"
  "explanation", ""
  "default", "Becke"
  "values", "(Becke|SSWeight)"


* grid/pruning

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "grid/pruning"
  "explanation", ""
  "default", "true"
  "values", "(yes|no)"


* TEI-integral-driven

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "TEI-integral-driven"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* J_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "J_engine"
  "explanation", ""
  "default", "RI_J"
  "values", "(conventional|RI_J|CD)"


* K_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "K_engine"
  "explanation", ""
  "default", "conventional"
  "values", "(conventional|RI_K|CD)"


* XC_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "XC_engine"
  "explanation", ""
  "default", "grid"
  "values", "(grid|gridfree|gridfree_CD)"


* CDAM_tau

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "CDAM_tau"
  "explanation", ""
  "default", "1.0E-10"
  "values", "(real)"


* CD_epsilon

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "CD_epsilon"
  "explanation", ""
  "default", "1.0E-4"
  "values", "(real)"


* CD/intermediate_file_format

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "CD/intermediate_file_format"
  "explanation", ""
  "default", "array_mmap"
  "values", "(array_mmap|array|mmap)"


* gridfree/orthogonalize_method

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "gridfree/orthogonalize_method"
  "explanation", ""
  "default", "canonical"
  "values", "(canonical|lowdin)"


* gridfree/dedicated_basis

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "gridfree/dedicated_basis"
  "explanation", ""
  "default", "false"
  "values", "(yes|no)"


* gridfree/save_v_eigval

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "gridfree/save_v_eigval"
  "explanation", ""
  "default", "true"
  "values", "(yes|no)"


* scf-memory-saving

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "scf-memory-saving"
  "explanation", ""
  "default", "no"
  "values", "yes/no"


* geometry/cartesian/input

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "geometry/cartesian/input"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* geometry/cartesian/unit

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "geometry/cartesian/unit"
  "explanation", ""
  "default", "a.u."
  "values", "au/a.u./angstrom"


* basis-set/orbital

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis-set/orbital"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* basis-set/density-auxiliary

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis-set/density-auxiliary"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* basis-set/exchange-auxiliary

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis-set/exchange-auxiliary"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* basis-set/gridfree

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis-set/gridfree"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* xc_density_threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "xc_density_threshold"
  "explanation", ""
  "default", "1.0E-16"
  "values", "real"


* geometry

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "geometry"
  "explanation", ""
  "default", "cartesian"
  "values", "cartesian/file"


* coordinates

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "coordinates"
  "explanation", ""
  "default", ""
  "values", ""


* basis_set

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis_set"
  "explanation", ""
  "default", ""
  "values", ""


* basis_set_j

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis_set_j"
  "explanation", ""
  "default", ""
  "values", ""


* basis_set_xc

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis_set_xc"
  "explanation", ""
  "default", ""
  "values", ""


* basis_set_gridfree

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "basis_set_gridfree"
  "explanation", ""
  "default", ""
  "values", ""


* independent-orbital-number

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "independent-orbital-number"
  "explanation", ""
  "default", "0"
  "values", "(integer > 0)"


* xc-potential/gxalpha/alpha-value

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "xc-potential/gxalpha/alpha-value"
  "explanation", ""
  "default", "0.7"
  "values", "real"


* myu-nyu-zero

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "myu-nyu-zero"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* guess/nsp-ppq

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/nsp-ppq"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* guess/sp-ppq

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/sp-ppq"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* guess/trans-angle-threshold

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/trans-angle-threshold"
  "explanation", ""
  "default", "1.0 1.5 20 30"
  "values", "nil/stored"


* guess/make-myu-nyu

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/make-myu-nyu"
  "explanation", ""
  "default", "meth0"
  "values", "meth0/meth1/meth2/meth3/meth4"


* guess/vct-normalize

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/vct-normalize"
  "explanation", ""
  "default", "ON OFF OFF"
  "values", ""


* guess/part-normalize

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/part-normalize"
  "explanation", ""
  "default", "nil"
  "values", "nil/stored"


* guess/user-vector

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "guess/user-vector"
  "explanation", ""
  "default", "nil"
  "values", "(nil|stored)"


* num_of_atoms

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_atoms"
  "explanation", ""
  "default", "0"
  "values", ""


* num_of_dummy_atoms

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_dummy_atoms"
  "explanation", ""
  "default", "0"
  "values", ""


* num_of_AOs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_AOs"
  "explanation", ""
  "default", "0"
  "values", ""


* num_of_MOs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_MOs"
  "explanation", ""
  "default", "0"
  "values", ""


* num_of_auxCDs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_auxCDs"
  "explanation", ""
  "default", "0"
  "values", ""


* num_of_auxXCs

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_auxXCs"
  "explanation", ""
  "default", "0"
  "values", ""


* TE

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "TE"
  "explanation", ""
  "default", ""
  "values", ""


* debug/file_warning

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/file_warning"
  "explanation", ""
  "default", "yes"
  "values", "(yes|no)"


* debug/save_J

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/save_J"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* debug/save_K

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/save_K"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* debug/save_Fxc_pure

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/save_Fxc_pure"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* debug/save_forces

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/save_forces"
  "explanation", ""
  "default", "no"
  "values", "(yes|no)"


* cutoff_distribution

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "cutoff_distribution"
  "explanation", ""
  "default", ""
  "values", ""


* length_scale_parameter

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "length_scale_parameter"
  "explanation", ""
  "default", "1"
  "values", ""


* control

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "control"
  "explanation", ""
  "default", ""
  "values", ""


* debug/eri/exact_J

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/eri/exact_J"
  "explanation", ""
  "default", ""
  "values", ""


* cutoff_density

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "cutoff_density"
  "explanation", ""
  "default", ""
  "values", ""


* cutoff_threshold_primitive

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "cutoff_threshold_primitive"
  "explanation", ""
  "default", "1e-12"
  "values", ""


* debug/eri/exact_K

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/eri/exact_K"
  "explanation", ""
  "default", ""
  "values", ""


* new_engine

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "new_engine"
  "explanation", ""
  "default", "true"
  "values", ""


* debug/eri/output_K

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/eri/output_K"
  "explanation", ""
  "default", ""
  "values", ""


* debug/eri/output_J

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/eri/output_J"
  "explanation", ""
  "default", ""
  "values", ""


* debug/DfXMatrix/save-eigval

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/DfXMatrix/save-eigval"
  "explanation", ""
  "default", "false"
  "values", ""


* debug/DfXMatrix/save-mat

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/DfXMatrix/save-mat"
  "explanation", ""
  "default", "false"
  "values", ""


* debug/DfXMatrix/check-X

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "debug/DfXMatrix/check-X"
  "explanation", ""
  "default", "false"
  "values", ""


* num_of_iterations

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "num_of_iterations"
  "explanation", ""
  "default", "0"
  "values", ""


* stat

.. csv-table::
  :widths: 20,80
  :stub-columns: 1

  "parameter", "stat"
  "explanation", ""
  "default", ""
  "values", ""



