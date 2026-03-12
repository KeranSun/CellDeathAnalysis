#' @keywords internal
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
  
  # 获取包版本
  version <- as.character(utils::packageVersion("CellDeathAnalysis"))
  
  # 创建启动消息
  msg <- paste0(
    "\n",
    "=======================================================================\n",
    "                                                                       \n",
    "     ___        _  _  ___                 _    _                       \n",
    "    / __| ___  | || ||   \\  ___  __ _ _ _| |_ | |_                     \n",
    "   | (__ / -_) | || || |) |/ -_)/ _` |  _|   \\|  _|                    \n",
    "    \\___|\\___| |_||_||___/ \\___|\\__,_|\\__|_||_|\\__|                    \n",
    "                                                                       \n",
    "         _                _           _                                \n",
    "        /_\\  _ _   __ _ | | _  _ ___(_) ___                            \n",
    "       / _ \\| ' \\ / _` || || || (_-<| |(_-<                            \n",
    "      /_/ \\_\\_||_|\\__,_||_| \\_, /__/|_|/__/                            \n",
    "                            |__/                                       \n",
    "=======================================================================\n",
    "\n",
    "  CellDeathAnalysis v", version, "\n",
    "  Comprehensive Analysis of Cell Death Pathways\n",
    "\n",
    "  Author: Keran Sun\n",
    "  Email:  s1214844197@163.com\n",
    "\n",
    "  Features: 14 cell death types | Scoring | Visualization\n",
    "            Survival | Enrichment | Single-cell | ML | Shiny App\n",
    "\n",
    "  Quick Start:\n",
    "    list_death_pathways()              # View all pathways\n",
    "    data(example_expr)                 # Load example data\n",
    "    launch_death_app()                 # Launch Shiny app\n",
    "\n",
    "  Documentation: ?CellDeathAnalysis\n",
    "\n",
    "=======================================================================\n"
  )
  
  packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
  # 可以在这里添加需要在加载时执行的代码
  invisible()
}
