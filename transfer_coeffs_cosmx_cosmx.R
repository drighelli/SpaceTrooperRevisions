## ============================================================
## 0. Libraries and output directory
## ============================================================
package_path <- "/Users/inzirio/My Drive/works/coding/SpaceTrooper"
devtools::load_all(package_path)
library(ggplot2)

out_dir <- "~/SpaceTrooper_qs_check1/cosmx_cosmx_transfer_check"
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
set.seed(1998)


## ============================================================
## 1. Define common CosMx model formula
## ============================================================
## Since both datasets are CosMx, we can use the CosMx-specific formula,
## including the border-adjusted aspect ratio term and control ratio.

cosmx_formula <- paste0(
    "~(",
    "log2SignalDensity + ",
    "Area_um + ",
    "I(abs(log2AspectRatio) * as.numeric(dist_border < 50)) + ",
    "log2Ctrl_total_ratio",
    ")^2"
)


## ============================================================
## 2. Read and preprocess CosMx Breast dataset
## ============================================================

cosmx_breast_path <- "~/Downloads/CosMx_data/DBKero/CosMx_Breast/CosMx_data_Case2/"

specosm_breast <- readCosmxSPE(cosmx_breast_path)

specosm_breast <- spatialPerCellQC(specosm_breast)


## ============================================================
## 3. Read and preprocess CosMx Pancreas dataset
## ============================================================

cosmx_pancreas_path <- "/Users/inzirio/Downloads/CosMx_data/Pancreas-CosMx-WTx"

specosm_pancreas <- readCosmxSPE(cosmx_pancreas_path)

specosm_pancreas <- spatialPerCellQC(specosm_pancreas)


## ============================================================
## 4. Train native QS models on both datasets
## ============================================================

specosm_breast <- computeQCScore(
    spe=specosm_breast,
    modelFormula=cosmx_formula,
    verbose=TRUE
)

specosm_pancreas <- computeQCScore(
    spe=specosm_pancreas,
    modelFormula=cosmx_formula,
    verbose=TRUE
)


## ============================================================
## 5. Store native QS scores and trained models
## ============================================================

specosm_breast$QC_score_breast_native <- specosm_breast$QC_score
specosm_pancreas$QC_score_pancreas_native <- specosm_pancreas$QC_score

breast_model <- metadata(specosm_breast)$QCScore_model
pancreas_model <- metadata(specosm_pancreas)$QCScore_model


## ============================================================
## 6. Apply transferred models
## ============================================================
## Apply Pancreas-trained model to Breast.
## Apply Breast-trained model to Pancreas.

specosm_breast <- applyQCScoreModel(
    spe=specosm_breast,
    qcModel=pancreas_model,
    scoreName="QC_score_pancreas_model"
)
saveRDS(specosm_breast, file.path(out_dir, "specosm_breast_with_pancreas_model.rds"))

specosm_pancreas <- applyQCScoreModel(
    spe=specosm_pancreas,
    qcModel=breast_model,
    scoreName="QC_score_breast_model"
)
saveRDS(specosm_pancreas, file.path(out_dir, "specosm_pancreas_with_breast_model.rds"))


## ============================================================
## 7. Basic checks
## ============================================================

summary(specosm_breast$QC_score_breast_native)
summary(specosm_breast$QC_score_pancreas_model)

summary(specosm_pancreas$QC_score_pancreas_native)
summary(specosm_pancreas$QC_score_breast_model)

range(specosm_breast$QC_score_breast_native, na.rm=TRUE)
range(specosm_breast$QC_score_pancreas_model, na.rm=TRUE)

range(specosm_pancreas$QC_score_pancreas_native, na.rm=TRUE)
range(specosm_pancreas$QC_score_breast_model, na.rm=TRUE)

sum(is.na(specosm_breast$QC_score_breast_native))
sum(is.na(specosm_breast$QC_score_pancreas_model))

sum(is.na(specosm_pancreas$QC_score_pancreas_native))
sum(is.na(specosm_pancreas$QC_score_breast_model))


## ============================================================
## 8. Correlations: native vs transferred
## ============================================================

cor_breast_spearman <- cor(
    specosm_breast$QC_score_breast_native,
    specosm_breast$QC_score_pancreas_model,
    use="complete.obs",
    method="spearman"
)

cor_breast_pearson <- cor(
    specosm_breast$QC_score_breast_native,
    specosm_breast$QC_score_pancreas_model,
    use="complete.obs",
    method="pearson"
)

cor_pancreas_spearman <- cor(
    specosm_pancreas$QC_score_pancreas_native,
    specosm_pancreas$QC_score_breast_model,
    use="complete.obs",
    method="spearman"
)

cor_pancreas_pearson <- cor(
    specosm_pancreas$QC_score_pancreas_native,
    specosm_pancreas$QC_score_breast_model,
    use="complete.obs",
    method="pearson"
)

cor_breast_spearman
cor_breast_pearson

cor_pancreas_spearman
cor_pancreas_pearson


## ============================================================
## 9. Plot: Breast native QS vs Pancreas-trained QS
## ============================================================
source("functions.R")
res_breast <- plot_qc_score_comparison(
    x=specosm_breast$QC_score_breast_native,
    y=specosm_breast$QC_score_pancreas_model,
    x_label="Breast native QS",
    y_label="Pancreas-trained QS",
    title="CosMx Breast: native QS vs Pancreas-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "cosmx_breast_native_vs_pancreas_model_A4_landscape.pdf"
    )
)

res_breast$plot
res_breast$correlation
#res_breast$quadrant_table

## ============================================================
## 10. Plot: Pancreas native QS vs Breast-trained QS
## ============================================================

res_pancreas <- plot_qc_score_comparison(
    x=specosm_pancreas$QC_score_pancreas_native,
    y=specosm_pancreas$QC_score_breast_model,
    x_label="Pancreas native QS",
    y_label="Breast-trained QS",
    title="CosMx Pancreas: native QS vs Breast-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "cosmx_pancreas_native_vs_breast_model_A4_landscape.pdf"
    )
)

res_pancreas$plot
res_pancreas$correlation
#res_pancreas$quadrant_table


## ============================================================
## 11. Compare model formulas and coefficients
## ============================================================

metadata(specosm_breast)$QCScore_model$model_formula
metadata(specosm_pancreas)$QCScore_model$model_formula

metadata(specosm_breast)$QCScore_model$model_matrix_colnames
metadata(specosm_pancreas)$QCScore_model$model_matrix_colnames

metadata(specosm_breast)$QCScore_model$coefficients_table
metadata(specosm_pancreas)$QCScore_model$coefficients_table


## ============================================================
## 12. Save summary table
## ============================================================

cosmx_transfer_summary <- data.frame(
    comparison=c(
        "Breast native vs Pancreas-trained",
        "Pancreas native vs Breast-trained"
    ),
    spearman=c(
        cor_breast_spearman,
        cor_pancreas_spearman
    ),
    pearson=c(
        cor_breast_pearson,
        cor_pancreas_pearson
    ),
    stringsAsFactors=FALSE
)

cosmx_transfer_summary

write.csv(
    cosmx_transfer_summary,
    file=file.path(out_dir, "cosmx_breast_pancreas_transfer_summary.csv"),
    row.names=FALSE
)

plot_qc_score_coldata <- function(x, y, spe, color_col,
    x_label="x", y_label="y",
    title="QC score comparison",
    threshold=0.75,
    cor_method="spearman",
    output_file=NULL) {

    stopifnot(length(x) == length(y))
    stopifnot(ncol(spe) == length(x))
    stopifnot(color_col %in% colnames(colData(spe)))

    df <- data.frame(
        x=x,
        y=y,
        color_var=as.data.frame(colData(spe))[[color_col]]
    )

    cor_val <- cor(
        df$x,
        df$y,
        use="complete.obs",
        method=cor_method
    )

    p <- ggplot(df, aes(x=x, y=y, color=color_var)) +
        geom_point(alpha=0.35, size=0.6) +
        geom_abline(
            slope=1,
            intercept=0,
            linetype="dashed",
            color="red"
        ) +
        geom_hline(yintercept=threshold) +
        geom_vline(xintercept=threshold) +
        labs(
            x=x_label,
            y=y_label,
            title=title,
            color=color_col
        ) +
        theme_minimal() +
        annotate(
            "text",
            x=quantile(df$x, 0.99, na.rm=TRUE),
            y=quantile(df$y, 0.01, na.rm=TRUE),
            label=paste0(
                cor_method,
                " r = ",
                formatC(cor_val, digits=3, format="f")
            ),
            hjust=1,
            vjust=0,
            color="black"
        )

    if (!is.null(output_file)) {
        ggsave(
            filename=output_file,
            plot=p,
            device="pdf",
            width=11.69,
            height=8.27,
            units="in"
        )
    }

    return(list(
        plot=p,
        correlation=cor_val,
        data=df
    ))
}

res_breast_signal <- plot_qc_score_coldata(
    x=specosm_breast$QC_score_breast_native,
    y=specosm_breast$QC_score_pancreas_model,
    spe=specosm_breast,
    color_col="log2SignalDensity",
    x_label="Breast native QS",
    y_label="Pancreas-trained QS",
    title="CosMx Breast: native QS vs Pancreas-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(out_dir, "cosmx_breast_native_vs_pancreas_model_colored_by_log2SignalDensity.pdf")
)

res_breast_signal$plot

res_breast_area <- plot_qc_score_coldata(
    x=specosm_breast$QC_score_breast_native,
    y=specosm_breast$QC_score_pancreas_model,
    spe=specosm_breast,
    color_col="Area_um",
    x_label="Breast native QS",
    y_label="Pancreas-trained QS",
    title="CosMx Breast: native QS vs Pancreas-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(out_dir, "cosmx_breast_native_vs_pancreas_model_colored_by_Area_um.pdf")
)

res_breast_area$plot

res_breast_ctrl <- plot_qc_score_coldata(
    x=specosm_breast$QC_score_breast_native,
    y=specosm_breast$QC_score_pancreas_model,
    spe=specosm_breast,
    color_col="log2Ctrl_total_ratio",
    x_label="Breast native QS",
    y_label="Pancreas-trained QS",
    title="CosMx Breast: native QS vs Pancreas-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(out_dir, "cosmx_breast_native_vs_pancreas_model_colored_by_log2Ctrl_total_ratio.pdf")
)

res_breast_ctrl$plot

res_breast_aspect <- plot_qc_score_coldata(
    x=specosm_breast$QC_score_breast_native,
    y=specosm_breast$QC_score_pancreas_model,
    spe=specosm_breast,
    color_col="log2AspectRatio",
    x_label="Breast native QS",
    y_label="Pancreas-trained QS",
    title="CosMx Breast: native QS vs Pancreas-trained QS",
    threshold=0.75,
    cor_method="spearman",
    output_file=file.path(out_dir, "cosmx_breast_native_vs_pancreas_model_colored_by_log2AspectRatio.pdf")
)

res_breast_aspect$plot

g1 <- SpaceTrooper::plotMetricHist(spe = specosm_breast, metric="log2Ctrl_total_ratio")
g2 <- SpaceTrooper::plotMetricHist(spe = specosm_pancreas, metric="log2Ctrl_total_ratio")
ggsave(
    filename=file.path(out_dir, "cosmx_breast_log2Ctrl_total_ratio_histogram.pdf"),
    plot=g1,
    device="pdf",
    width=11.69,
    height=8.27,
    units="in"
)
ggsave(
    filename=file.path(out_dir, "cosmx_pancreas_log2Ctrl_total_ratio_histogram.pdf"),
    plot=g2,
    device="pdf",
    width=11.69,
    height=8.27,
    units="in"
)
############################################################
specosm_mb1 <- readCosmxSPE("/Users/inzirio/Downloads/CosMx_data/CosMx1k_MouseBrain1")
specosm_mb2 <- readCosmxSPE("/Users/inzirio/Downloads/CosMx_data/CosMx1k_MouseBrain2")


specosm_mb1 <- spatialPerCellQC(specosm_mb1)
specosm_mb2 <- spatialPerCellQC(specosm_mb2)

cosmx_formula <- paste0(
    "~(",
    "log2SignalDensity + ",
    "Area_um + ",
    "I(abs(log2AspectRatio) * as.numeric(dist_border < 50)) + ",
    "log2Ctrl_total_ratio",
    ")^2"
)

specosm_mb1 <- computeQCScore(
    spe=specosm_mb1,
    modelFormula=cosmx_formula,
    verbose=TRUE
)

specosm_mb2 <- computeQCScore(
    spe=specosm_mb2,
    modelFormula=cosmx_formula,
    verbose=TRUE
)

specosm_mb1$QC_score_mb1_native <- specosm_mb1$QC_score
specosm_mb2$QC_score_mb2_native <- specosm_mb2$QC_score

mb1_model <- metadata(specosm_mb1)$QCScore_model
mb2_model <- metadata(specosm_mb2)$QCScore_model

saveRDS(specosm_mb1, file.path(out_dir, "specosm_mb1_with_native_qs.rds"))
saveRDS(specosm_mb2, file.path(out_dir, "specosm_mb2_with_native_qs.rds"))
specosm_mb1 <- applyQCScoreModel(
    spe=specosm_mb1,
    qcModel=mb2_model,
    scoreName="QC_score_mb2_model"
)

specosm_mb2 <- applyQCScoreModel(
    spe=specosm_mb2,
    qcModel=mb1_model,
    scoreName="QC_score_mb1_model"
)

cor_mb1_spearman <- cor(
    specosm_mb1$QC_score_mb1_native,
    specosm_mb1$QC_score_mb2_model,
    use="complete.obs",
    method="spearman"
)

cor_mb1_pearson <- cor(
    specosm_mb1$QC_score_mb1_native,
    specosm_mb1$QC_score_mb2_model,
    use="complete.obs",
    method="pearson"
)

cor_mb2_spearman <- cor(
    specosm_mb2$QC_score_mb2_native,
    specosm_mb2$QC_score_mb1_model,
    use="complete.obs",
    method="spearman"
)

cor_mb2_pearson <- cor(
    specosm_mb2$QC_score_mb2_native,
    specosm_mb2$QC_score_mb1_model,
    use="complete.obs",
    method="pearson"
)

res_mb1 <- plot_qc_score_comparison(
    x=specosm_mb1$QC_score_mb1_native,
    y=specosm_mb1$QC_score_mb2_model,
    x_label="MouseBrain1 native QS",
    y_label="MouseBrain2-trained QS",
    title="CosMx MouseBrain1: native QS vs MouseBrain2-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "cosmx_mousebrain1_native_vs_mousebrain2_model_A4_landscape.pdf"
    )
)

res_mb2 <- plot_qc_score_comparison(
    x=specosm_mb2$QC_score_mb2_native,
    y=specosm_mb2$QC_score_mb1_model,
    x_label="MouseBrain2 native QS",
    y_label="MouseBrain1-trained QS",
    title="CosMx MouseBrain2: native QS vs MouseBrain1-trained QS",
    threshold=NULL, #change to 0.75 if you want to visualize quadrants
    cor_method="spearman",
    output_file=file.path(
        out_dir,
        "cosmx_mousebrain2_native_vs_mousebrain1_model_A4_landscape.pdf"
    )
)

res_mb1$plot
res_mb2$plot

color_vars <- c(
    "log2SignalDensity",
    "Area_um",
    "log2AspectRatio",
    "log2Ctrl_total_ratio"
)

mb1_colored_plots <- lapply(color_vars, function(v) {
    plot_qc_score_coldata(
        x=specosm_mb1$QC_score_mb1_native,
        y=specosm_mb1$QC_score_mb2_model,
        spe=specosm_mb1,
        color_col=v,
        x_label="MouseBrain1 native QS",
        y_label="MouseBrain2-trained QS",
        title=paste0("MouseBrain1 native vs MouseBrain2-trained QS: ", v),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("cosmx_mousebrain1_native_vs_mousebrain2_model_", v, ".pdf")
        )
    )
})
names(mb1_colored_plots) <- color_vars
mb1_colored_plots
mb2_colored_plots <- lapply(color_vars, function(v) {
    plot_qc_score_coldata(
        x=specosm_mb2$QC_score_mb2_native,
        y=specosm_mb2$QC_score_mb1_model,
        spe=specosm_mb2,
        color_col=v,
        x_label="MouseBrain2 native QS",
        y_label="MouseBrain1-trained QS",
        title=paste0("MouseBrain2 native vs MouseBrain1-trained QS: ", v),
        threshold=0.75,
        cor_method="spearman",
        output_file=file.path(
            out_dir,
            paste0("cosmx_mousebrain2_native_vs_mousebrain1_model_", v, ".pdf")
        )
    )
})
names(mb2_colored_plots) <- color_vars
mb2_colored_plots

cosmx_mousebrain_transfer_summary <- data.frame(
    comparison=c(
        "MouseBrain1 native vs MouseBrain2-trained",
        "MouseBrain2 native vs MouseBrain1-trained"
    ),
    spearman=c(cor_mb1_spearman, cor_mb2_spearman),
    pearson=c(cor_mb1_pearson, cor_mb2_pearson),
    stringsAsFactors=FALSE
)

cosmx_mousebrain_transfer_summary

write.csv(
    cosmx_mousebrain_transfer_summary,
    file=file.path(out_dir, "cosmx_mousebrain1_mousebrain2_transfer_summary.csv"),
    row.names=FALSE
)

metadata(specosm_mb1)$QCScore_model$model_matrix_colnames
metadata(specosm_mb2)$QCScore_model$model_matrix_colnames

####################
