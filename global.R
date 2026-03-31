## ============================================================
## GLOBAL SETUP
## ============================================================

## Required packages — install manually if missing, e.g. install.packages(c("shiny", ...))
local({
  required_packages <- c(
    "shiny", "later", "DT", "htmlwidgets",
    "shinyjs", "shinyWidgets",
    "pedtools", "kinship2", "pedmut", "forrel",
    "pedFamilias",
    "tibble", "dplyr", "purrr", "stringr"
  )
  missing_pkgs <- required_packages[!vapply(
    required_packages,
    function(p) requireNamespace(p, quietly = TRUE),
    logical(1)
  )]
  if (length(missing_pkgs)) {
    stop(
      "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
      ". Install with:\n  install.packages(c(",
      paste0('"', missing_pkgs, '"', collapse = ", "),
      "), dependencies = TRUE)",
      call. = FALSE
    )
  }
  for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
  }

  # If user/session enables {conflicted}, make app intent explicit for known overlaps.
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(dplyr::filter, dplyr::bind_rows)
  }
})

## Internal helpers live on the search path, not in .GlobalEnv
KINSEARCH <- "KinshipAssembly"
if (KINSEARCH %in% search()) {
  detach(KINSEARCH, character.only = TRUE)
}
.kinship_fn_env <- new.env(parent = globalenv())
sys.source("logic/functions.R", envir = .kinship_fn_env, keep.source = FALSE)
attach(.kinship_fn_env, pos = 2L, name = KINSEARCH)
rm(.kinship_fn_env, KINSEARCH)

## Shiny options
options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB max upload per request

## When deployed (e.g. Shiny Server / Connect), avoid leaking stack traces to the browser
if (nzchar(Sys.getenv("SHINY_PORT", "")) || nzchar(Sys.getenv("RSTUDIO_PRODUCT", ""))) {
  options(shiny.sanitize.errors = TRUE)
}

## Application version
VERSION <- "1.0"

## Global configuration
## Logging: INFO on when running locally, off when `SHINY_PORT` is set (typical deployment).
## Override with KINSHIP_LOG_INFO / KINSHIP_LOG_DEBUG = true, false, 1, or 0 (case-insensitive).
.env_flag <- function(name, default) {
  v <- trimws(Sys.getenv(name, ""))
  if (!nzchar(v)) return(default)
  v <- tolower(v)
  if (v %in% c("1", "true", "yes", "t")) return(TRUE)
  if (v %in% c("0", "false", "no", "f")) return(FALSE)
  default
}
.local_run <- !nzchar(Sys.getenv("SHINY_PORT", ""))
CONFIG <- list(
  mp_id = "Missing person",
  ## Bridge label only: EXTRA_1 -> this -> mp_id before mergePed (never use EXTRA_ as merge key)
  merge_mp_transient_label = "__KINSHIP_MP__",
  mut_model = "none",
  mut_rate = 0.002,
  exclude_patterns = "-(RM|RP)",
  info = if (nzchar(trimws(Sys.getenv("KINSHIP_LOG_INFO", "")))) {
    .env_flag("KINSHIP_LOG_INFO", .local_run)
  } else {
    .local_run
  },
  debug = if (nzchar(trimws(Sys.getenv("KINSHIP_LOG_DEBUG", "")))) {
    .env_flag("KINSHIP_LOG_DEBUG", FALSE)
  } else {
    FALSE
  },
  default_markers = NULL
)
rm(.local_run, .env_flag)

## UI configuration (titles, metadata)
UI_CONFIG <- list(
  app_title = "KinshipAssembly",
  app_subtitle = "Genetic Comparison between MPI/DVI and POI Components Using STR Markers",
  version = VERSION,
  author = "sbiagini0"
)

## Custom CSS (see www/style.css). Inline script: non-detail modals use max-content width;
## the detail modal (kinship-detail-modal) uses a fixed width in CSS — here we only run
## column reflow for the table and clear obsolete inline widths.
## shiny::icon() uses Font Awesome bundled with Shiny (no extra CSS required).
custom_css <- tags$head(
  tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
  tags$script(HTML(
"
(function() {
  function horizPadBorder(el) {
    if (!el) return 0;
    var s = window.getComputedStyle(el);
    return parseFloat(s.paddingLeft) + parseFloat(s.paddingRight) +
      parseFloat(s.borderLeftWidth) + parseFloat(s.borderRightWidth);
  }
  function kinshipModalIsResetConfirm(modal) {
    if (!modal || !modal.querySelector) return false;
    if (modal.querySelector('.kinship-reset-confirm-marker')) return true;
    var dlg = modal.querySelector('.modal-dialog');
    return !!(dlg && dlg.classList.contains('kinship-reset-confirm-modal'));
  }
  function kinshipClearModalInlineWidths(modal) {
    if (!modal || !modal.querySelector) return;
    var dlg = modal.querySelector('.modal-dialog');
    var mc = modal.querySelector('.modal-content');
    var mb = modal.querySelector('.modal-body');
    var wrap = modal.querySelector('.dataTables_wrapper');
    var tbl = wrap && wrap.querySelector ? wrap.querySelector('table') : null;
    [dlg, mc, mb].forEach(function(el) {
      if (!el) return;
      el.style.removeProperty('width');
      el.style.removeProperty('min-width');
      el.style.removeProperty('max-width');
    });
    if (wrap) {
      wrap.style.removeProperty('width');
      wrap.style.removeProperty('max-width');
    }
    if (tbl) tbl.style.removeProperty('width');
  }
  function kinshipPlotBlockMinWidth(modal) {
    var block = modal.querySelector('.kinship-modal-pedigree-block');
    if (!block) return 0;
    var plotOut = block.querySelector('.shiny-plot-output');
    if (!plotOut) return 0;
    var w = Math.max(plotOut.offsetWidth, plotOut.scrollWidth);
    var img = plotOut.querySelector('img');
    if (img && img.offsetWidth) w = Math.max(w, img.offsetWidth);
    if (img && img.complete && img.naturalWidth) w = Math.max(w, img.naturalWidth);
    return w;
  }
  function kinshipFitWideDataTableModal(modal) {
    if (!modal || !modal.querySelector) return;
    var dlg = modal.querySelector('.modal-dialog');
    var mc = modal.querySelector('.modal-content');
    var mb = modal.querySelector('.modal-body');
    var wrap = modal.querySelector('.dataTables_wrapper');
    var tbl = wrap && wrap.querySelector ? wrap.querySelector('table') : null;
    /* Detail modal: fixed width in CSS only. */
    if (dlg && dlg.classList.contains('kinship-detail-modal')) {
      kinshipClearModalInlineWidths(modal);
      return;
    }
    /* Reset confirm: never run wide-table sizing (marker + class on .modal-dialog). */
    if (kinshipModalIsResetConfirm(modal)) {
      kinshipClearModalInlineWidths(modal);
      return;
    }
    var hasTable = !!(wrap && tbl);
    [dlg, mc, mb].forEach(function(el) {
      if (!el) return;
      el.style.removeProperty('width');
      el.style.removeProperty('min-width');
      el.style.removeProperty('max-width');
    });
    if (wrap) {
      wrap.style.removeProperty('width');
      wrap.style.removeProperty('max-width');
    }
    if (tbl) tbl.style.removeProperty('width');
    void modal.offsetHeight;
    var wPlot = kinshipPlotBlockMinWidth(modal);
    var sw = mc ? mc.scrollWidth : 0;
    if (hasTable) {
      wrap.style.setProperty('width', 'max-content', 'important');
      wrap.style.setProperty('box-sizing', 'border-box', 'important');
      tbl.style.setProperty('width', 'auto', 'important');
      void modal.offsetHeight;
      sw = Math.max(sw, mc.scrollWidth);
      if (sw < 200) {
        var wTbl = Math.max(
          wrap.scrollWidth,
          tbl.scrollWidth,
          tbl.offsetWidth
        );
        sw = Math.max(sw, wTbl + horizPadBorder(wrap) + horizPadBorder(mb) + horizPadBorder(mc));
      }
    }
    /* Reserve space for modal-body padding + pedigree inset (see .kinship-modal-pedigree-block CSS) */
    var mbHorizPad = 0;
    if (mb) {
      var csmb = window.getComputedStyle(mb);
      mbHorizPad = parseFloat(csmb.paddingLeft) + parseFloat(csmb.paddingRight);
    }
    var pedigreeInset = 28;
    sw = Math.max(sw, wPlot + 48 + pedigreeInset);
    var inner = Math.max(sw + 32, wPlot > 0 ? wPlot + 64 + pedigreeInset : 0, 400);
    /* Extra slack so refit passes do not shrink the dialog flush to the plot (margin “disappears”) */
    var safety = 48;
    var w = Math.ceil(Math.min(inner + mbHorizPad + safety, window.innerWidth * 0.96));
    if (w < 320) w = Math.ceil(Math.min(1000, window.innerWidth * 0.96));
    if (dlg) {
      dlg.style.setProperty('width', w + 'px', 'important');
      dlg.style.setProperty('max-width', '96vw', 'important');
    }
    if (mc) {
      mc.style.setProperty('width', w + 'px', 'important');
      mc.style.setProperty('max-width', '96vw', 'important');
      mc.style.setProperty('box-sizing', 'border-box', 'important');
      mc.style.setProperty('overflow', 'visible', 'important');
    }
    if (mb) {
      mb.style.setProperty('width', '100%', 'important');
      mb.style.setProperty('max-width', '100%', 'important');
      mb.style.setProperty('box-sizing', 'border-box', 'important');
      mb.style.setProperty('overflow-x', 'visible', 'important');
    }
    if (wrap && hasTable) {
      wrap.style.setProperty('width', '100%', 'important');
      wrap.style.setProperty('max-width', '100%', 'important');
      wrap.style.setProperty('box-sizing', 'border-box', 'important');
    }
  }
  window.kinshipFitWideDataTableModal = kinshipFitWideDataTableModal;
  /* Detail modal: on first click the table may not exist yet or width is 0 during the fade;
     retry until DataTables is ready and after several layout ticks. */
  function kinshipEqualizeDetailModalTable(modal) {
    if (!window.jQuery) return false;
    /* Do not require .modal.show: drawCallback may run during the fade, before .show */
    if (!modal) modal = document.querySelector('#shiny-modal');
    if (!modal) return false;
    var wrap = modal.querySelector('#lr_details_modal_table');
    if (!wrap) return false;
    var tbl = wrap.querySelector('table');
    if (!tbl || !jQuery.fn.DataTable || !jQuery.fn.DataTable.isDataTable(tbl)) return false;
    var api = jQuery(tbl).DataTable();
    var n = api.columns().count();
    if (n < 1) return false;
    var pct = (100 / n) + '%';
    var css = { width: pct, minWidth: 0, maxWidth: 'none', boxSizing: 'border-box' };
    jQuery(api.table().node()).css({ width: '100%', tableLayout: 'fixed' });
    for (var i = 0; i < n; i++) {
      jQuery(api.column(i).header()).css(css);
      jQuery(api.column(i).nodes()).css(css);
    }
    return true;
  }
  window.kinshipEqualizeDetailModalTable = kinshipEqualizeDetailModalTable;
  function kinshipTryEqualizeDetailModal(modal, retriesLeft) {
    if (!modal) return;
    if (kinshipEqualizeDetailModalTable(modal)) {
      jQuery(window).trigger('resize');
      return;
    }
    if (retriesLeft > 0) {
      setTimeout(function() { kinshipTryEqualizeDetailModal(modal, retriesLeft - 1); }, 50);
    }
  }
  function kinshipDetailModalBootstrap(modal) {
    if (!modal || !modal.querySelector) return;
    var dlg = modal.querySelector('.modal-dialog.kinship-detail-modal');
    if (!dlg) return;
    kinshipTryEqualizeDetailModal(modal, 60);
    [0, 50, 100, 200, 400, 600].forEach(function(ms) {
      setTimeout(function() {
        kinshipEqualizeDetailModalTable(modal);
        jQuery(window).trigger('resize');
      }, ms);
    });
  }
  function kinshipRefitOpenModal() {
    var modal = document.querySelector('#shiny-modal.modal.show');
    if (!modal) return;
    if (kinshipModalIsResetConfirm(modal)) return;
    kinshipFitWideDataTableModal(modal);
  }
  var kinshipModalResizeObserver = null;
  var kinshipResizeDebounce = null;
  function kinshipScheduleModalRefit() {
    kinshipRefitOpenModal();
    setTimeout(kinshipRefitOpenModal, 0);
    setTimeout(kinshipRefitOpenModal, 50);
    setTimeout(kinshipRefitOpenModal, 150);
    setTimeout(kinshipRefitOpenModal, 350);
    setTimeout(kinshipRefitOpenModal, 700);
    setTimeout(kinshipRefitOpenModal, 1200);
  }
  if (window.jQuery) {
    jQuery(document).on('shown.bs.modal', function(e) {
      var modal = e.target;
      if (kinshipModalIsResetConfirm(modal)) {
        kinshipClearModalInlineWidths(modal);
      } else {
        kinshipFitWideDataTableModal(modal);
        kinshipScheduleModalRefit();
      }
      var dlgDetail = modal.querySelector('.modal-dialog.kinship-detail-modal');
      if (dlgDetail) {
        kinshipDetailModalBootstrap(modal);
      }
      var mc = modal.querySelector('.modal-content');
      var dlgEl = modal.querySelector('.modal-dialog');
      /* Do not observe detail / reset modals: refit on resize is for wide-table dialogs only */
      var skipRO = dlgEl && (dlgEl.classList.contains('kinship-detail-modal') ||
        kinshipModalIsResetConfirm(modal));
      if (mc && typeof ResizeObserver !== 'undefined' && !skipRO) {
        if (kinshipModalResizeObserver) kinshipModalResizeObserver.disconnect();
        kinshipModalResizeObserver = new ResizeObserver(function() {
          if (kinshipResizeDebounce) clearTimeout(kinshipResizeDebounce);
          kinshipResizeDebounce = setTimeout(function() {
            kinshipResizeDebounce = null;
            kinshipRefitOpenModal();
          }, 80);
        });
        kinshipModalResizeObserver.observe(mc);
      }
    });
    jQuery(document).on('hidden.bs.modal', function() {
      if (kinshipResizeDebounce) {
        clearTimeout(kinshipResizeDebounce);
        kinshipResizeDebounce = null;
      }
      if (kinshipModalResizeObserver) {
        kinshipModalResizeObserver.disconnect();
        kinshipModalResizeObserver = null;
      }
    });
  }
  if (window.jQuery) {
    jQuery(document).on('shiny:value', function(ev) {
      var t = ev.target;
      if (!t || !t.id) return;
      if (t.id === 'hypothesis_plot_modal' || t.id === 'lr_details_modal_table') {
        kinshipScheduleModalRefit();
      }
      if (t.id === 'lr_details_modal_table') {
        var modalOpen = document.querySelector('#shiny-modal');
        if (modalOpen) kinshipDetailModalBootstrap(modalOpen);
      }
      if (t.id === 'hypothesis_plot_modal') {
        setTimeout(function() { jQuery(window).trigger('resize'); }, 0);
        setTimeout(function() { jQuery(window).trigger('resize'); }, 120);
      }
    });
  }
  if (window.jQuery) {
    jQuery(document).on('shown.bs.modal', function() {
      var imgs = document.querySelectorAll('#shiny-modal.modal.show .kinship-modal-pedigree-block img');
      imgs.forEach(function(img) {
        if (!img.complete) {
          img.addEventListener('load', function onload() {
            img.removeEventListener('load', onload);
            kinshipRefitOpenModal();
            setTimeout(kinshipRefitOpenModal, 50);
          });
        }
      });
    });
  }
})();
"
  ))
)
