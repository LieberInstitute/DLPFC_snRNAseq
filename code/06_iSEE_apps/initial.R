initial <- list()

################################################################################
# Settings for Reduced dimension plot 1
################################################################################

initial[["ReducedDimensionPlot1"]] <- new("ReducedDimensionPlot",
    Type = "UMAP", XAxis = 1L, YAxis = 2L,
    FacetRowByColData = "SAMPLE_ID", FacetColumnByColData = "SAMPLE_ID",
    ColorByColumnData = "cellType_layer", ColorByFeatureNameAssay = "logcounts",
    ColorBySampleNameColor = "#FF0000", ShapeByColumnData = "SAMPLE_ID",
    SizeByColumnData = "age", FacetRowBy = "None", FacetColumnBy = "None",
    ColorBy = "Column data", ColorByDefaultColor = "#000000",
    ColorByFeatureName = "MIR1302-2HG", ColorByFeatureSource = "---",
    ColorByFeatureDynamicSource = FALSE, ColorBySampleName = "1_AAACCCAAGTTCTCTT-1",
    ColorBySampleSource = "---", ColorBySampleDynamicSource = FALSE,
    ShapeBy = "None", SizeBy = "None", SelectionAlpha = 0.1,
    ZoomData = numeric(0), BrushData = list(), VisualBoxOpen = FALSE,
    VisualChoices = "Color", ContourAdd = FALSE, ContourColor = "#0000FF",
    PointSize = 1, PointAlpha = 1, Downsample = FALSE, DownsampleResolution = 200,
    CustomLabels = FALSE, CustomLabelsText = "1_AAACCCAAGTTCTCTT-1",
    FontSize = 1, LegendPointSize = 1, LegendPosition = "Bottom",
    HoverInfo = TRUE, LabelCenters = FALSE, LabelCentersBy = "SAMPLE_ID",
    LabelCentersColor = "#000000", VersionInfo = list(iSEE = structure(list(
        c(2L, 8L, 0L)
    ), class = c("package_version", "numeric_version"))), PanelId = c(ReducedDimensionPlot = 1L), PanelHeight = 600L,
    PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list()
)

################################################################################
# Settings for Complex heatmap 1
################################################################################

initial[["ComplexHeatmapPlot1"]] <- new("ComplexHeatmapPlot",
    Assay = "logcounts", CustomRows = TRUE,
    CustomRowsText = "CLDN5 \nGFAP \nGAD1\nSLC17A7 \nTMEM119\nOLIG2\nSNAP25\nMBP\nPDGFRA",
    ClusterRows = FALSE, ClusterRowsDistance = "spearman", ClusterRowsMethod = "ward.D2",
    DataBoxOpen = FALSE, VisualChoices = "Annotations", ColumnData = c(
        "cellType_layer",
        "Sample"
    ), RowData = character(0), CustomBounds = FALSE,
    LowerBound = NA_real_, UpperBound = NA_real_, AssayCenterRows = FALSE,
    AssayScaleRows = FALSE, DivergentColormap = "purple < black < yellow",
    ShowDimNames = "Rows", LegendPosition = "Right", LegendDirection = "Vertical",
    VisualBoxOpen = FALSE, NamesRowFontSize = 10, NamesColumnFontSize = 10,
    ShowColumnSelection = TRUE, OrderColumnSelection = TRUE,
    VersionInfo = list(iSEE = structure(list(c(2L, 8L, 0L)), class = c(
        "package_version",
        "numeric_version"
    ))), PanelId = c(ComplexHeatmapPlot = 1L),
    PanelHeight = 600L, PanelWidth = 8L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)

################################################################################
# Settings for Row data table 1
################################################################################

initial[["RowDataTable1"]] <- new("RowDataTable",
    Selected = "SNAP25", Search = "SNA", SearchColumns = c(
        "",
        "", "", "", "", "", "", "", ""
    ), HiddenColumns = character(0),
    VersionInfo = list(iSEE = structure(list(c(2L, 8L, 0L)), class = c(
        "package_version",
        "numeric_version"
    ))), PanelId = c(RowDataTable = 1L), PanelHeight = 600L,
    PanelWidth = 4L, SelectionBoxOpen = FALSE, RowSelectionSource = "---",
    ColumnSelectionSource = "---", DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE,
    ColumnSelectionDynamicSource = FALSE, RowSelectionRestrict = FALSE,
    ColumnSelectionRestrict = FALSE, SelectionHistory = list()
)

################################################################################
# Settings for Feature assay plot 1
################################################################################

initial[["FeatureAssayPlot1"]] <- new("FeatureAssayPlot",
    Assay = "logcounts", XAxis = "Column data",
    XAxisColumnData = "cellType_layer", XAxisFeatureName = "MIR1302-2HG",
    XAxisFeatureSource = "---", XAxisFeatureDynamicSource = FALSE,
    YAxisFeatureName = "SNAP25", YAxisFeatureSource = "RowDataTable1",
    YAxisFeatureDynamicSource = TRUE, FacetRowByColData = "SAMPLE_ID",
    FacetColumnByColData = "SAMPLE_ID", ColorByColumnData = "cellType_layer",
    ColorByFeatureNameAssay = "logcounts", ColorBySampleNameColor = "#FF0000",
    ShapeByColumnData = "SAMPLE_ID", SizeByColumnData = "age",
    FacetRowBy = "None", FacetColumnBy = "None", ColorBy = "Column data",
    ColorByDefaultColor = "#000000", ColorByFeatureName = "MIR1302-2HG",
    ColorByFeatureSource = "---", ColorByFeatureDynamicSource = FALSE,
    ColorBySampleName = "1_AAACCCAAGTTCTCTT-1", ColorBySampleSource = "---",
    ColorBySampleDynamicSource = FALSE, ShapeBy = "None", SizeBy = "None",
    SelectionAlpha = 0.1, ZoomData = numeric(0), BrushData = list(),
    VisualBoxOpen = FALSE, VisualChoices = "Color", ContourAdd = FALSE,
    ContourColor = "#0000FF", PointSize = 1, PointAlpha = 1,
    Downsample = FALSE, DownsampleResolution = 200, CustomLabels = FALSE,
    CustomLabelsText = "1_AAACCCAAGTTCTCTT-1", FontSize = 1,
    LegendPointSize = 1, LegendPosition = "Bottom", HoverInfo = TRUE,
    LabelCenters = FALSE, LabelCentersBy = "SAMPLE_ID", LabelCentersColor = "#000000",
    VersionInfo = list(iSEE = structure(list(c(2L, 8L, 0L)), class = c(
        "package_version",
        "numeric_version"
    ))), PanelId = c(FeatureAssayPlot = 1L),
    PanelHeight = 600L, PanelWidth = 8L, SelectionBoxOpen = FALSE,
    RowSelectionSource = "---", ColumnSelectionSource = "---",
    DataBoxOpen = FALSE, RowSelectionDynamicSource = FALSE, ColumnSelectionDynamicSource = FALSE,
    RowSelectionRestrict = FALSE, ColumnSelectionRestrict = FALSE,
    SelectionHistory = list()
)
