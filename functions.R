
# To be decided if this should be moved to a separate file or included in the main package
plot_qc_score_comparison <- function(x, y, x_label="x", y_label="y",
    title="QC score comparison", threshold=NULL,
    cor_method="spearman", output_file=NULL) {

    stopifnot(length(x) == length(y))

    df <- data.frame(x=x, y=y)

    # If `threshold` is NULL or non-finite (e.g. Inf), disable quadrant logic
    use_threshold <- !(is.null(threshold) || !is.finite(threshold))

    if (use_threshold) {
        df$quadrant <- with(df, ifelse(
            x < threshold & y < threshold, "low / low",
            ifelse(x >= threshold & y < threshold, "high x / low y",
            ifelse(x < threshold & y >= threshold, "low x / high y",
                "high / high"))
        ))
    } else {
        df$quadrant <- "all"
    }

    cor_val <- cor(df$x, df$y, use="complete.obs", method=cor_method)

    # Build the ggplot object depending on whether the threshold is used.
    if (use_threshold) {
        p <- ggplot(df, aes(x=x, y=y, color=quadrant)) +
            geom_point(alpha=0.35, size=0.6) +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            geom_hline(yintercept=threshold) +
            geom_vline(xintercept=threshold) +
            scale_color_manual(values=c(
                "low / low"="grey60",
                "high x / low y"="orange",
                "low x / high y"="dodgerblue",
                "high / high"="black"
            )) +
            labs(x=x_label, y=y_label, title=title, color="Quadrant") +
            theme_minimal()
    } else {
        p <- ggplot(df, aes(x=x, y=y)) +
            geom_point(alpha=0.35, size=0.6) +
            geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
            labs(x=x_label, y=y_label, title=title) +
            theme_minimal()
    }

    p <- p + annotate(
        "text",
        x=quantile(df$x, 0.99, na.rm=TRUE),
        y=quantile(df$y, 0.01, na.rm=TRUE),
        label=paste0(cor_method, " r = ", formatC(cor_val, digits=3, format="f")),
        hjust=1,
        vjust=0,
        color="black"
    )

    if (!is.null(output_file)) {
        ggsave(output_file, p, device="pdf",
            width=11.69, height=8.27, units="in")
    }

    return(list(
        plot=p,
        correlation=cor_val,
        quadrant_table=table(df$quadrant)
    ))
}
